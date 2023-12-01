#' @name import_coverage_scenario
#' @title import_coverage_scenario
#' @description Imports a CSV file with vaccination coverage values for the scenario
#' @param modelpath Path to montagu files
#' @param country country code
#' @param scenario Unique string that identifies the coverage scenario name
#' @param filter0 logical for whether 0 vaccination years should be filtered out (default = FALSE)
#' @param redownload whether to redownload the file
#' @importFrom magrittr %>%
#' @return Dataframe with vaccination coverage for a single scenario and country. Years without vaccination are excluded from the returned dataframe. Countries without vaccination in any year return a null dataframe.
#' @export
#' @include retrieve_montagu_coverage.R
import_coverage_scenario <- function(modelpath, country, scenario, num_doses = NULL, filter0 = FALSE, redownload = TRUE){
  
  ##calam added to ensure csv montagu coverage for each vaccination scenario is pulled from the montagu folder for the 202310gavi-4 touchstone 
  runname <- config$runname
  if (runname == "202310gavi-4" & scenario == "campaign-default"){
    num_doses <- config$vacc$ndoses
  }
  if (runname == "202310gavi-4" & scenario == "campaign-default" & num_doses == "one"){
    scenario <- "ocv1-default"
  } else if (runname == "202310gavi-4" & scenario == "campaign-default" & num_doses == "two"){
    scenario <- "ocv1-ocv2-default"
  }
  ##end addition
  
  #First check, then retrieve
  CoverageFiles <- list.files(modelpath, pattern = "^coverage_")
  if (length(CoverageFiles) >= 2 & redownload == FALSE){ #there should be 2 coverage data files, so the default should be 2
    message(paste0("The coverage data files have been under the directory: ", modelpath, '. No new download was made. '))
  } else{
    retrieve_montagu_coverage(modelpath)
  }
  
  #Start importing
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
#' @description Imports a CSV template file for the central burden estimates static columns indicating the burden estimates that need to be generated for each scenario
#' @param modelpath Path to montagu files
#' @param country country code
#' @param redownload whether to redownload the file
#' @return dataframe for central burden template for one country
#' @export 
#' @include retrieve_montagu_centralburden_template.R
import_centralburden_template <- function(modelpath, country, redownload = TRUE){
  
  #First check, then retrieve
  CentralBurdenTempFiles <- list.files(modelpath, pattern = "^central-burden")
  if (length(CentralBurdenTempFiles) > 0 & redownload == FALSE){
    message(paste0("The central burden template files have been under the directory: ", modelpath, '. No new download was made. '))
  } else{
    retrieve_montagu_centralburden_template(modelpath)
  }
  
  #Start importing
  fname <- list.files(modelpath, pattern = "^central-burden")
  if (length(fname)>1){
    stop(paste("More than 1 central burden template was found in", modelpath))
  }
  rc <- readr::read_csv(paste0(modelpath, "/", fname)) %>%
    dplyr::filter(country == !!country) %>%
    dplyr::distinct(disease, country, country_name, year, age)
  return(rc)
}


#' @name import_country_population
#' @title import_country_population
#' @description Get total country population from standardized source in Montagu. [may be used in future for case-based targeting]
#' @param modelpath path to montagu files
#' @param country country code
#' @param redownload whether to redownload the file
#' @importFrom magrittr %>%
#' @return 
#' @export 
#' @include retrieve_montagu_population.R
import_country_population <- function(modelpath, country, redownload = TRUE){
  
  #First check, then retrieve
  TotPopFiles <- list.files(modelpath, pattern = "tot_pop_both.csv$")
  if (length(TotPopFiles) >= 1 & redownload == FALSE){ #there should be 1 file
    message(paste0("The total population files have been under the directory: ", modelpath, '. No new download was made. '))
  } else{
    retrieve_montagu_population(modelpath)
  }
  
  #Start importing
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

  country_pop <- import_country_population(modelpath, country, redownload = FALSE) #this is important
  if (year > max(country_pop$year)){
    message(paste("Population data not available for", year, "- Using", max(country_pop$year), "data instead"))
    year <- max(country_pop$year)
  }
  year_pop <- country_pop[which(country_pop$year == as.numeric(year)),]$pop_model

  return(year_pop)
}

#' @name import_country_agePop
#' @title import_country_agePop
#' @description Get age-specific (1-year) country population from standardized source in Montagu. [may be used in future for case-based targeting]
#' @param modelpath path to montagu files
#' @param country country code
#' @param redownload whether to redownload the file
#' @importFrom magrittr %>%
#' @return 
#' @export 
#' @include retrieve_montagu_agePop.R
import_country_agePop <- function(modelpath, country, redownload = TRUE){
  
  #First check, then retrieve
  AgePopFiles <- list.files(modelpath, pattern = "int_pop_both.csv$")
  if (length(AgePopFiles) >= 1 & redownload == FALSE){ #there should be 1 file
    message(paste0("The age-specific population files have been under the directory: ", modelpath, '. No new download was made. '))
  } else{
    retrieve_montagu_agePop(modelpath) 
  }
  
  #Start importing
  pop_fn <- list.files(modelpath, pattern = "int_pop_both.csv$")
  if (length(pop_fn)>1){
    stop(paste("More than 1 int_pop_both demographic file was found in", modelpath))
  }
  message(paste0("Loading ", modelpath, "/", pop_fn))
  agepop <- readr::read_csv(paste0(modelpath, "/", pop_fn)) %>%
    dplyr::filter(country_code == !!country) %>%
    dplyr::rename(GID_0 = country_code, pop_age = value) %>%
    dplyr::select(GID_0, year, age_from, age_to, pop_age)

  totpop <- dplyr::group_by(agepop, GID_0, year) %>%
    dplyr::summarise(pop_tot = sum(pop_age))

  rc <- dplyr::left_join(agepop, totpop, by = c("GID_0", "year")) %>%
    dplyr::mutate(prop_age = pop_age/pop_tot) %>%
    dplyr::rename(country = GID_0)

  return(rc)
}


#' @name import_country_lifeExpectancy
#' @title import_country_lifeExpectancy
#' @description Get total life expectancy at birth from standardized source in Montagu. [may be used in future for case-based targeting]
#' @param modelpath path to montagu files
#' @param country country code
#' @param redownload whether to redownload the file
#' @importFrom magrittr %>%
#' @return 
#' @export 
#' @include retrieve_montagu_lifeExpectancy.R
import_country_lifeExpectancy <- function(modelpath, country, redownload = TRUE){
  
  #First check, then retrieve
  LifeExpFiles <- list.files(modelpath, pattern = "lx0_both.csv$")
  if (length(LifeExpFiles) >= 1 & redownload == FALSE){ #there should be 1 file
    message(paste0("The life expectancy data files have been under the directory: ", modelpath, '. No new download was made. '))
  } else{
    retrieve_montagu_lifeExpectancy(modelpath)
  }
  
  #Start importing
  lx0_fn <- list.files(modelpath, pattern = "lx0_both.csv$")
  if (length(lx0_fn)>1){
    stop(paste("More than 1 lx0_both (life expectancy) file was found in", modelpath))
  }
  message(paste0("Loading ", modelpath, "/", lx0_fn))
  rc <- readr::read_csv(paste0(modelpath, "/", lx0_fn)) %>%
    dplyr::filter(country_code == !!country) %>%
    dplyr::select(country_code, year, value) %>%
    dplyr::rename(country = country_code, lx0 = value)

  return(rc)
}

#' @name import_country_lifeExpectancy_1yr
#' @title import_country_lifeExpectancy_1yr
#' @description Get total life expectancy at birth for a single year from standardized source in Montagu. [may be used in future for case-based targeting]
#' @param modelpath path to montagu files
#' @param country country code
#' @param year year for which life expectancy should be returned
#' @importFrom magrittr %>%
#' @return
#' @export
import_country_lifeExpectancy_1yr <- function(modelpath, country, year){
  lx0_df <- import_country_lifeExpectancy(modelpath, country, redownload = FALSE) #important change
  rc <- dplyr::filter(lx0_df, year == !!year) %>%
    dplyr::select(lx0) %>%
    unlist %>% unname

  return(rc)
}


#' @name import_disability_weight
#' @title import_disability_weight
#' @description Returns mean disability weight for severe diarrheal disease from IHME Global Burden of Disease 2016 and 2019 (same estimate)
#' @return numeric value for disability weight
#' @export 
import_disability_weight <- function(){
  ## 95% CI 0.164 to 0.348 
  ## See input_data/IHME_GBD_2019_DISABILITY_WEIGHTS_Y2020M010D15
  return(0.247)
}

#' @name import_templateFilename_prefix
#' @title import_templateFilename_prefix
#' @description Return string from filename template
#' @param type "stochastic", "central", "parameter"
#' @param modelpath path to Montagu files
#' @param redownload whether to redownload the file
#' @return template filename prefix
#' @export 
#' @include retrieve_montagu_stochasticburden_template.R
import_templateFilename_prefix <- function(type, modelpath, redownload = FALSE){
  
  #First check, then retrieve
  StocBurdenTempFiles <- list.files(modelpath, pattern = "stochastic-burden-template")
  if (length(StocBurdenTempFiles) > 0 & redownload == FALSE){
    message(paste0("The stochastic burden template file has been under the directory: ", modelpath, '. No new download was made. '))
  } else{
    retrieve_montagu_stochasticburden_template(modelpath)
  }
  
  #The file is already available
  if (type == "stochastic"){
    fn <- list.files(modelpath, pattern = "stochastic-burden-template")[1]
    rc <- stringr::str_remove(fn, " standard template.csv")
  } else if (type == "central"){
    fn <- list.files(modelpath, pattern = "central-burden-template")[1]
    rc <- stringr::str_remove(fn, " standard template.csv")
  } else if (type == "parameter"){
    fn <- list.files(modelpath, pattern = "stochastic-burden-template")[1]
    rc <- stringr::str_replace(stringr::str_remove(fn, " standard template.csv"), "stochastic-burden", "stochastic-params")
  } else{
    warning(paste("This", type, "filename is not supported."))
    rc <- NULL
  }

  return(rc)

}


#' @name import_country_proportion_under5
#' @title import_country_proportion_under5
#' @description Get proportion of country population aged 1-4 from standardized source in Montagu
#'
#' @param modelpath path to montagu files
#' @param country country code
#' @param year year (year in which vaccination campaign takes place)
#' @param redownload whether to redownload the file
#'
#' @importFrom magrittr %>%
#' @return the proportion of the population aged <5 per year for the specified country and year
#' @export 
#' @include retrieve_montagu_agePop.R
import_country_proportion_under5 <- function(modelpath, country, year, redownload = TRUE){
  
  #First check, then retrieve
  AgePopFiles <- list.files(modelpath, pattern = "int_pop_both.csv$")
  if (length(AgePopFiles) >= 1 & redownload == FALSE){ #there should be 1 file
    message(paste0("The age-specific population files have been under the directory: ", modelpath, '. No new download was made. '))
  } else{
    retrieve_montagu_agePop(modelpath) 
  }
  
  #Start importing
  pop_fn <- list.files(modelpath, pattern = "int_pop_both.csv$")
  if (length(pop_fn)>1){
    stop(paste("More than 1 int_pop_both demographic file was found in", modelpath))
  }
  message(paste0("Loading ", modelpath, "/", pop_fn))
  agepop <- readr::read_csv(paste0(modelpath, "/", pop_fn)) %>%
    dplyr::filter(country_code == !!country) %>%
    dplyr::filter(!age_from %in% c(0,100)) %>%   ##this filters out the people aged 0
    dplyr::rename(GID_0 = country_code, pop_age = value) %>%
    dplyr::select(GID_0, year, age_from, age_to, pop_age)
  
  totpop <- dplyr::group_by(agepop, GID_0, year) %>%
    dplyr::summarise(pop_tot = sum(pop_age))
  
  rc <- dplyr::left_join(agepop, totpop, by = c("GID_0", "year")) %>%
    dplyr::mutate(prop_age = pop_age/pop_tot) %>%
    dplyr::rename(country = GID_0) 
    
    ##this is a major point differentiating this function from import_country_agePop
    ##get proportion of the population that is under 5 years old for each year
    under5_allyears <- rc %>%
    dplyr::filter(age_from %in% c(1,2,3,4)) %>%
    dplyr::group_by(year) %>%
    summarise(prop_under5 = sum(prop_age)) 
  
  under5 <- under5_allyears[under5_allyears$year == year,]$prop_under5
  
  return(under5)
}

