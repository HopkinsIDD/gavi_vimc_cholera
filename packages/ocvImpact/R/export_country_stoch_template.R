#' @name export_country_stoch_template
#' @title export_country_stoch_template
#' @description  Export a stochastic burden template file for a single country
#' @param modelpath path to montagu files
#' @param country country code
#' @param scenario Unique string that identifies the coverage scenario name
#' @param rawoutpath path to raw model output files
#' @param outpath path to final model output files
#' @importFrom magrittr %>%
#' @return 
#' @export
export_country_stoch_template <- function(
  modelpath,
  country, 
  scenario,
  rawoutpath,
  outpath
  ){

  ## import templates
  cb_template <- import_centralburden_template(modelpath, country)

  ## calculation inputs
  disab_wt <- import_disability_weight()
  infect_dur <- generate_infectionDuration()
  cfr <- generate_cfr(country)
  aoi <- generate_aoi(country) ## average age of infection
  lifeExpect_df <- import_country_lifeExpectancy(modelpath, country)

  if(!all(unique(cb_template$year) %in% lifeExpect_df$year)){
    missing_yrs <- unique(cb_template$year)[which(!unique(cb_template$year) %in% lifeExpect_df$year)]
    adhoc_row <- dplyr::filter(lifeExpect_df, year == missing_yrs-1)
    adhoc_row$year <- missing_yrs
    lifeExpect_df <- dplyr::bind_rows(lifeExpect_df, adhoc_row) %>%
      dplyr::arrange(country, year)
    warning(paste("Replacing missing life expectancy data from", missing_yrs, "with", missing_yrs-1, "data, respectively."))
  }

  ## import already-generated model outputs
  ec_out_fn <- paste0(rawoutpath, "/", scenario, "/", country, "_ec.csv")
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
  pop_age_df <- import_country_agePop(modelpath, country)
  stoch <- dplyr::left_join(expCases_age, pop_age_df, by = c("country", "year")) %>%
    dplyr::mutate(
      cases = round(cases_tot*prop_age, 0),
      deaths = round(deaths_tot*prop_age, 0),
      yll = round(yll_tot*prop_age, 0),
      yld = round(yld_tot*prop_age, 0),
      dalys = round(daly_tot*prop_age, 0)
      ) %>%
    dplyr::rename(age = age_from, cohort_size = pop_age) %>%
    dplyr::full_join(cb_template, by = c("year", "country", "age")) %>%
    dplyr::select(disease, run_id, year, age, country, country_name, cohort_size, cases, deaths, dalys) %>%
    dplyr::arrange(run_id, age, year)

  out_fn <- paste0(outpath, "/", country, "_", scenario, "_stoch.csv")
  message(paste("Writing", out_fn))
  readr::write_csv(stoch, out_fn)

  ## get stochastic parameters
  params <- dplyr::distinct(expCases, run_id, country, incid_rate) %>%
    dplyr::mutate(aoi = aoi, cfr = cfr, infect_dur = infect_dur) %>%
    dplyr::select(run_id, country, aoi, incid_rate, cfr, infect_dur) 
  par_fn <- paste0(outpath, "/", country, "_", scenario, "_pars.csv")
  message(paste("Writing", par_fn))
  readr::write_csv(params, par_fn)

  return(stoch)
  
}
