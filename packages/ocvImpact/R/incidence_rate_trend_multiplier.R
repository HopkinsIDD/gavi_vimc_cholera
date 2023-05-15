#' @name incidence_rate_trend_multiplier
#' @title incidence_rate_trend_multiplier
#' @description Create and return a function that can take in the year to give a prediction of the incidence rate in that country
#' @param datapath path to input data 
#' @param modelpath path to montagu files
#' @param country country code
#' @param datathreshold the threshold for the minimum number of years of data in the WHO dataset to give a valid incidence rate estimate (the default is 5)
#' @param use_country_incid_trend whether we wanna use the country-level who data for a certain country. It overrides others. 
#' @return a function that can take in years and give a incidence rate prediction for a given country
#' @export
#' @include utils.R utils_montagu.R
incidence_rate_trend_multiplier <- function(datapath,
                                            modelpath,
                                            country,
                                            datathreshold = 5, 
                                            use_country_incid_trend){
  
  ### get the WHO file first
  who_dir_file_name <- paste0(datapath, '/incidence/who_case_repo_source.csv')
  if (file.exists(who_dir_file_name)){
    casedf <- readr::read_csv(who_dir_file_name)
  }else{
    stop('The WHO country-level annual incidence dataset does not exist under the correct directory. Please check. ')
  }
  
  ### it seems necessary to load the dplyr package here -- 11/26
  library(dplyr)

  ### get the population and join two datasets
  SplittedString = strsplit(modelpath, '/')[[1]]
  touchstone = SplittedString[length(SplittedString)]
  tot_pop <- montagu::montagu_demographic_data("tot_pop", touchstone) %>% 
    select(country_code, year, value) %>% 
    mutate(pop_model = value) %>% 
    select( - value)
  
  casedf_pop <- casedf %>% 
    mutate(country_code = ISO3) %>% 
    select(-ISO3)
  casedf_pop <- left_join(casedf_pop, tot_pop, by = c("country_code", "year")) %>% 
    mutate(ic_rate = sCh/pop_model)
  print('The check for three selects in the incidence function has been passed. ')

  ### clean the dataset and run two regression models -- country and general 
  casedf_pop <- casedf_pop[!is.na(casedf_pop$sCh), ] %>% 
    mutate(offset = pop_model * 1)
  trend.country <- glm(sCh ~ year 
                            + as.factor(country_code) 
                            + year * as.factor(country_code)
                            + offset(log(offset)), 
                            data = casedf_pop, family = quasipoisson)
  trend.general <- glm(sCh ~ year 
                           + offset(log(offset)), 
                           data = casedf_pop, family = quasipoisson)
  
  ### get the specific country estimate and make a function
  if (country %in% casedf_pop$country_code 
      & length(unique(casedf_pop[casedf_pop$country_code == country, ]$year)) >= as.numeric(datathreshold)
      & use_country_incid_trend){
    ## If we're satisfied with the WHO data for this country
    incidence_rate_trend_function <- function(year, country_code = country) {
      ### The country might be the reference country
      if (country == sort(unique(casedf_pop$country_code))[1]){
        return(as.numeric(
          (exp(trend.country$coefficients["(Intercept)"] + 
              as.numeric(year) * (trend.country$coefficients["year"]))) /
            (exp(trend.country$coefficients["(Intercept)"] + 
                   2014 * (trend.country$coefficients["year"])))
          ))
      }else{
        return(as.numeric(
          (exp(trend.country$coefficients["(Intercept)"] + 
              trend.country$coefficients[paste0('as.factor(country_code)', country)] +
              as.numeric(year) * (trend.country$coefficients["year"] + 
                                  trend.country$coefficients[paste0('year:as.factor(country_code)', country)]))) / 
            (exp(trend.country$coefficients["(Intercept)"] + 
                   trend.country$coefficients[paste0('as.factor(country_code)', country)] +
                   2014 * (trend.country$coefficients["year"] + 
                                         trend.country$coefficients[paste0('year:as.factor(country_code)', country)])))
          ))
      }
        
    }
  }else{
    ## If we are not satisfied with the WHO data for the country or just don't wanna use anyways
    warning('WHO annual incidence dataset does not have any data or enough data for this country, the global mean estimate will be used instead. ')
    incidence_rate_trend_function <- function(year, country_code = country) {
      return(as.numeric(
        (exp(trend.general$coefficients["(Intercept)"] + 
            as.numeric(year) * (trend.general$coefficients["year"]))) / 
          (exp(trend.general$coefficients["(Intercept)"] + 
                 2014 * (trend.general$coefficients["year"])))
      ))
      
    }
  }
  return(incidence_rate_trend_function)
  
}

