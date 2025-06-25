# 25 June 2025, Qulu Zheng & Elizabeth Lee
# This script is used to process the final output from the VIMC core models with careseeking, community death, and targeitng adjustments

# Functions ----------------------------------------------------------------

#' @name adjust_cases
#' @title adjust_cases
#' @description adjust the cases in the final model output file
#' @param df the final model output file
#' @param use_stochastic whether to adjust the model output in a stochastic way
#' @param multiplier the single multiplier for cases, the default is 1 (un-adjusted)
#' @return a dataframe with adjusted cases (same format as the final model output file)
#' @export
adjust_cases <- function(df, use_stochastic = F, multiplier = 1) {
  multiplier <- if (!use_stochastic) multiplier else stop("stochastic scaling for cases is not implemented yet.")
  
  df$cases <- df$cases * multiplier
  
  return(df)
}

#' @name adjust_deaths
#' @title adjust_deaths
#' @description adjust the deaths in the final model output file
#' @param df the final model output file
#' @param use_stochastic whether to adjust the model output in a stochastic way
#' @param multiplier the single multiplier for deaths, the default is 1 (un-adjusted)
#' @return a dataframe with adjusted deaths (same format as the final model output file)
#' @export

adjust_deaths <- function(df, use_stochastic = F, multiplier = 1) {
  multiplier <- if (!use_stochastic) multiplier else stop("stochastic scaling for deaths is not implemented yet.")
  
  df$deaths <- df$deaths * multiplier
  
  return(df)
}

#' @name adjust_DALYs
#' @title adjust_DALYs
#' @description adjust the DALYs (YLLs + YLDs) in the final model output file
#' @param df the final model output file
#' @param use_stochastic whether to adjust the cases and deaths in the model output in a stochastic way
#' @param deaths_multiplier the single multiplier for deaths, the default is 1 (un-adjusted)
#' @param cases_multiplier the single multiplier for cases, the default is 1 (un-adjusted)
#' @return a dataframe with adjusted DALYs (same format as the final model output file)
#' @export
adjust_DALYs <- function(df, use_stochastic = F, cases_multiplier = 1, deaths_multiplier = 1) {
  cases_multiplier <- if (!use_stochastic) cases_multiplier else stop("stochastic scaling for cases is not implemented yet.")
  deaths_multiplier <- if (!use_stochastic) deaths_multiplier else stop("stochastic scaling for deaths is not implemented yet.")
  
  df$yld <- (df$dalys-df$yll) * cases_multiplier
  df$yll <- df$yll * deaths_multiplier
  df$dalys <- df$yld + df$yll
  
  return(df)
}

# Prepare paths and load data ---------------------------------------------------------------
filepath_novacc <- #<path to original novaccination stochastic model output>
filepath_vacc1 <- #<path to original ocv1 stochastic model output>
filepath_vacc2 <- #<path to original ocv1-ocv2 stochastic model output>
  
new_filepath_novacc <- #<path to new novaccination stochastic model output>
new_filepath_vacc1 <- #<path to new ocv1 stochastic model output>
new_filepath_vacc2 <- #<path to new ocv1-ocv2 stochastic model output>

# Set parameters ----------------------------------------------------
targeting_adjustment_on <- TRUE # logical indicating whether to include adjustment for targeting strategy

case_mult_novacc <- 1/.328 # multiplier for careseeking adjustment
death_mult_novacc <- 3.87 # multiplier for community death to facility death ratio
targeting_mult <- 1/1.7 # multiplier for changing from affected population to MAI-based targeting

case_mult_vacc <- ifelse(targeting_adjustment_on, case_mult_novacc * targeting_mult, case_mult_novacc)
death_mult_vacc <- ifelse(targeting_adjustment_on, death_mult_novacc * targeting_mult, death_mult_novacc)

# Adjust model outputs ----------------------------------------------------
novacc <- readr::read_csv(filepath_novacc)
vacc1 <- readr::read_csv(filepath_vacc1)
vacc2 <- readr::read_csv(filepath_vacc2)

novacc_adj <- novacc %>% 
  adjust_cases(., multiplier = case_mult_novacc) %>% 
  adjust_deaths(., multiplier  = death_mult_novacc) %>% 
  adjust_DALYs(., cases_multiplier = case_mult_novacc, deaths_multiplier = death_mult_novacc)

vacc1_adj <- vacc1 %>% 
  adjust_cases(., multiplier = case_mult_vacc) %>% 
  adjust_deaths(., multiplier  = death_mult_vacc) %>% 
  adjust_DALYs(., cases_multiplier = case_mult_vacc, deaths_multiplier = death_mult_vacc)

vacc2_adj <- vacc2 %>% 
  adjust_cases(., multiplier = case_mult_vacc) %>% 
  adjust_deaths(., multiplier  = death_mult_vacc) %>% 
  adjust_DALYs(., cases_multiplier = case_mult_vacc, deaths_multiplier = death_mult_vacc)

# Write adjusted model outputs ----------------------------------------------------
readr::write_csv(new_filepath_novacc)
readr::write_csv(new_filepath_vacc1)
readr::write_csv(new_filepath_vacc2)