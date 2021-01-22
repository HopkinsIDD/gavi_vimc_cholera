#' @name export_scenario_stoch_template
#' @title export_scenario_stoch_template
#' @description Save a stochastic burden template file for all countries in a scenario
#' @param modelpath path to montagu files
#' @param country country code
#' @param scenario Unique string that identifies the coverage scenario name
#' @param rawoutpath path to raw model output files
#' @return 
#' @export
export_country_stoch_template <- function(
  modelpath,
  country, 
  scenario,
  rawoutpath
  ){

  ## import templates
  cb_template <- import_centralburden_template(modelpath, country)

  ## calculation inputs
  disab_wt <- import_disability_weight()
  infect_dur <- generate_infectionDuration()
  cfr <- generate_cfr(country)
  aoi <- generate_aoi(country) ## average age of infection
  lifeExpect_df <- import_country_lifeExpectancy(modelpath, country)

  ## import already-generated model outputs
  ec_out_fn <- paste0(rawoutpath, "/", scenario, "/", country, "_ec.csv")
  expCases <- readr::read_csv(ec_out_fn) %>%
    dplyr::left_join(lifeExpect_df, by = c("country", "year")) %>%
    dplyr::mutate(
      deaths_tot = cfr*ec,
      aoi_cl = min(c(aoi, lx0)), ## set aoi to lower value between life expectancy and aoi table
      yll_tot = ed*(lx0-aoi_cl),
      yld_tot = ec*infect_dur*disab_wt,
      daly_tot = yll+yld
      ) %>%
    dplyr::rename(cases_tot = ec)
  
  ## distribute model outputs by proportion of the population
  pop_age <- import_country_agePop(modelpath, country)
  stoch <- dplyr::left_join(pop_age, expCases, by = c("country", "year")) %>%
    dplyr::mutate(
      cases = round(cases_tot*prop_age, 0),
      deaths = round(deaths_tot*prop_age, 0),
      yll = round(yll_tot*prop_age, 0),
      yld = round(yld_tot*prop_age, 0),
      daly = round(daly_tot*prop_age, 0)
      ) %>%
    dplyr::rename(age = age_from, cohort_size = agepop) %>%
    dplyr::select(disease, run_id, year, age, country, country_name, cohort_size, cases, deaths, dalys) %>%
    dplyr::arrange(run_id, age, year)

  return(stoch)
}