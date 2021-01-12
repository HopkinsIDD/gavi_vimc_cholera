#' @name import_coverage_scenario
#' @title import_coverage_scenario
#' @description Imports a CSV file with vaccination coverage values for the scenario
#' @param modelpath Path to montagu files
#' @param country country code
#' @param scenario Unique string that identifies the coverage scenario name
#' @param filter0 logical for whether 0 vaccination years should be filtered out (default = FALSE)
#' @importFrom magrittr %>%
#' @return Dataframe with vaccination coverage for a single scenario and country. Years without vaccination are excluded from the returned dataframe. Countries without vaccination in any year return a null dataframe.
#' @export
import_coverage_scenario <- function(modelpath, country, scenario, filter0 = FALSE){
  cov_fnames <- list.files(modelpath, pattern = "^coverage_")
  scn_fname <- paste0(modelpath, "/", grep(scenario, cov_fnames, ignore.case = TRUE, value = TRUE))

  if (length(scn_fname)>1){
    stop(paste("You have not specified a unique scenario in", modelpath, "with the following string:", scenario))
  }

  message(paste("Loading coverage scenario:", country, scenario))
  cov_dat <- readr::read_csv(scn_fname) %>%
    dplyr::filter(country_code == !!country) 

  if (nrow(cov_dat)==0){
    message("**** N.B. This is a no-vaccination scenario. Returning NULL.")
    cov_dat <- NULL
  } else{
    cov_dat <- dplyr::arrange(cov_dat, year) %>%
      dplyr::mutate(fvp = round(target*coverage, 0))
    if (filter0){
      cov_dat <- dplyr::filter(cov_dat, coverage != 0)
    }
  }

  return(cov_dat)
}

#' @name import_centralburden_template
#' @title import_centralburden_template
#' @description Imports a CSV template file for the central burden estimates indicating the burden estimates that need to be generated for each scenario
#' @param modelpath Path to montagu files
#' @param country country code
#' @return dataframe for central burden template for one country
#' @export 
import_centralburden_template <- function(modelpath, country){
  fname <- list.files(modelpath, pattern = "^central-burden")
  if (length(fname)>1){
    stop(paste("More than 1 central burden template was found in", modelpath))
  }
  readr::read_csv(paste0(modelpath, "/", fname)) %>%
    dplyr::filter(country == country)
}


#' @name import_country_population
#' @title import_country_population
#' @description Get total country population from standardized source in Montagu. [may be used in future for case-based targeting]
#' @param modelpath path to montagu files
#' @param country country code
#' @importFrom magrittr %>%
#' @return 
#' @export 
import_country_population <- function(modelpath, country){
  pop_fn <- list.files(modelpath, pattern = "tot_pop_both.csv$")
  if (length(pop_fn)>1){
    stop(paste("More than 1 tot_pop_both demographic file was found in", modelpath))
  }
  message(paste0("Loading ", modelpath, "/", pop_fn))
  rc <- readr::read_csv(paste0(modelpath, "/", pop_fn)) %>%
    dplyr::filter(country_code == !!country) %>%
    dplyr::rename(GID_0 = country_code, pop_model = value) %>%
    dplyr::select(GID_0, year, pop_model)

  return(rc)
}


#' @name import_country_population_1yr
#' @title import_country_population_1yr
#' @description Get total country population from standardized source in Montagu for a single year
#' @param modelpath path to montagu files
#' @param country country code
#' @param year population year
#' @return numeric value
#' @export
import_country_population_1yr <- function(modelpath, country, year){

  country_pop <- import_country_population(modelpath, country)
  year_pop <- country_pop[which(country_pop$year == as.numeric(year)),]$pop_model

  return(year_pop)
}
