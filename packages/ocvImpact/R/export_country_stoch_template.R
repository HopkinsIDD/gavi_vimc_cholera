#' @name export_country_stoch_template
#' @title export_country_stoch_template
#' @description  Export a stochastic burden template file for a single country
#' @param modelpath path to montagu files
#' @param country country code
#' @param scenario Unique string that identifies the coverage scenario name
#' @param rawoutpath path to raw model output files
#' @param outpath path to final model output files
#' @importFrom magrittr %>%
#' @return dataframe
#' @export
#' @include utils.R utils_montagu.R
export_country_stoch_template <- function(
  modelpath,
  country,
  scenario,
  rawoutpath,
  outpath
  ){

  ## avoid reading montagu files multiple times
  montagu_cache <- new.env()
  montagu_cache[["centralburden_template"]] <- import_centralburden_template(modelpath, country, montagu_cache, redownload = FALSE)
  montagu_cache[["country_agePop"]] <- import_country_agePop(modelpath, country, montagu_cache, redownload = FALSE)
  montagu_cache[["country_lifeExpectancy"]] <- import_country_lifeExpectancy(modelpath, country, montagu_cache, redownload = FALSE)


  ## import templates
  cb_template <- montagu_cache[["centralburden_template"]]

  ## calculation inputs
  disab_wt <- import_disability_weight()
  infect_dur <- generate_infectionDuration()
  cfr <- generate_cfr(country)
  aoi <- generate_aoi(country) ## average age of infection
  lifeExpect_df <- montagu_cache[["country_lifeExpectancy"]]

  if(!all(unique(cb_template$year) %in% lifeExpect_df$year)){
    missing_yrs <- unique(cb_template$year)[which(!unique(cb_template$year) %in% lifeExpect_df$year)]
    adhoc_row <- dplyr::filter(lifeExpect_df, year == missing_yrs-1)
    adhoc_row$year <- missing_yrs
    lifeExpect_df <- dplyr::bind_rows(lifeExpect_df, adhoc_row) %>%
      dplyr::arrange(country, year)
    warning(paste("Replacing missing life expectancy data from", missing_yrs, "with", missing_yrs-1, "data, respectively."))
  }

  ## import already-generated model outputs
  incidence_rate_trend <- as.logical(config$setting$incidence_rate_trend)
  outbreak_multiplier <- as.logical(config$setting$outbreak_multiplier)
  setting <- paste0('incid_trend_', incidence_rate_trend, '_outb_layer_',  outbreak_multiplier)
  ##calam added for the new one dose and two dose campaigns for the 202310gavi-4 touchstone
  runname <- config$runname
  if (runname == '202310gavi-4'){
    ndoses <- config$vacc$ndoses
    ec_out_fn <- paste0(rawoutpath, "/", scenario, "/", setting, "/", country,"_",ndoses,"_ec.csv")
  } else {
    ec_out_fn <- paste0(rawoutpath, "/", scenario, "/", setting, "/", country, "_ec.csv")
  }
  ##calam addition end
  expCases <- readr::read_csv(ec_out_fn)
  expCases_age <- dplyr::left_join(expCases, lifeExpect_df, by = c("country", "year")) %>%
    dplyr::mutate(
      ed = cfr*ec,
      aoi_cl = pmin(aoi, lx0), ## set aoi to lower value between life expectancy and aoi table
      yll_tot = ed*(lx0-aoi_cl),
      yld_tot = ec*infect_dur*disab_wt,
      daly_tot = yll_tot+yld_tot
      ) %>%
    dplyr::rename(cases_tot = ec, deaths_tot = ed)

  ## distribute model outputs by proportion of the population
  pop_age_df <- montagu_cache[["country_agePop"]]
  ##calam added to reflect change in central burden template for touchstone 202310gavi-4, which requires yll
  if (runname == '202310gavi-4'){
    stoch <- dplyr::left_join(expCases_age, pop_age_df, by = c("country", "year")) %>%
      dplyr::mutate(
        cases = round(cases_tot*prop_age, 4),
        deaths = round(deaths_tot*prop_age, 4),
        yll = round(yll_tot*prop_age, 4),
        yld = round(yld_tot*prop_age, 4),
        dalys = round(daly_tot*prop_age, 4)
      ) %>%
      dplyr::rename(age = age_from, cohort_size = pop_age) %>%
      dplyr::full_join(cb_template, by = c("year", "country", "age")) %>%
      dplyr::select(disease, run_id, year, age, country, country_name, cohort_size, cases, deaths, dalys, yll) %>%
      dplyr::arrange(run_id, age, year)
  } else {
    stoch <- dplyr::left_join(expCases_age, pop_age_df, by = c("country", "year")) %>%
      dplyr::mutate(
        cases = round(cases_tot*prop_age, 4),
        deaths = round(deaths_tot*prop_age, 4),
        yll = round(yll_tot*prop_age, 4),
        yld = round(yld_tot*prop_age, 4),
        dalys = round(daly_tot*prop_age, 4)
      ) %>%
      dplyr::rename(age = age_from, cohort_size = pop_age) %>%
      dplyr::full_join(cb_template, by = c("year", "country", "age")) %>%
      dplyr::select(disease, run_id, year, age, country, country_name, cohort_size, cases, deaths, dalys) %>%
      dplyr::arrange(run_id, age, year)
  }
  ##end addition

  ## include setting into the file name
  # incidence_rate_trend <- as.logical(config$setting$incidence_rate_trend)
  # outbreak_multiplier <- as.logical(config$setting$outbreak_multiplier)
  # setting <- paste0('incid_trend_', incidence_rate_trend, '_outb_layer_',  outbreak_multiplier)

  ##calam added added for the new one dose and two dose campaigns for the 202310gavi-4 touchstone
  if (runname == '202310gavi-4'){
    ndoses <- config$vacc$ndoses
    out_fn <- paste0(outpath, "/", country, "_", scenario, "_", ndoses, "_", setting, "_stoch.csv")
  } else {
    out_fn <- paste0(outpath, "/", country, "_", scenario, "_", setting, "_stoch.csv")
  }
  ##calam addition end
  message(paste("Writing", out_fn))
  readr::write_csv(stoch, out_fn)

  ## get stochastic parameters
  params <- expCases %>%
    dplyr::group_by(run_id, country) %>%
    dplyr::summarise(incid_rate = mean(incid_rate)) %>%
    dplyr::mutate(aoi = aoi, cfr = cfr, infect_dur = infect_dur) %>%
    dplyr::select(run_id, country, aoi, incid_rate, cfr, infect_dur)
  ##calam added added for the new one dose and two dose campaigns for the 202310gavi-4 touchstone
  if (runname == '202310gavi-4'){
    ndoses <- config$vacc$ndoses
    par_fn <- paste0(outpath, "/", country, "_", scenario, "_", ndoses, "_", setting, "_pars.csv")
  }else {
    par_fn <- paste0(outpath, "/", country, "_", scenario, "_", setting, "_pars.csv")
  }
  ##calam addition end
  message(paste("Writing", par_fn))
  readr::write_csv(params, par_fn)

  return(stoch)

}
