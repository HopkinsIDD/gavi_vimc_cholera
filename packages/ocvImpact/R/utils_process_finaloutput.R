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
#' @description adjust the DALYs(YLLs) in the final model output file
#' @param df the final model output file
#' @param deaths_use_stochastic whether to adjust the deaths in the model output in a stochastic way
#' @param cases_use_stochastic whether to adjust the cases in the model output in a stochastic way
#' @param deaths_multiplier the single multiplier for deaths, the default is 1 (un-adjusted)
#' @param cases_multiplier the single multiplier for cases, the default is 1 (un-adjusted)
#' @return a dataframe with adjusted DALYs(YLLs) (same format as the final model output file)
#' @export
adjust_DALYs <- function(df, cases_use_stochastic = F, deaths_use_stochastic = F, cases_multiplier = 1, deaths_multiplier = 1) {
  cases_multiplier <- if (!cases_use_stochastic) cases_multiplier else stop("stochastic scaling for deaths is not implemented yet.")
  deaths_multiplier <- if (!deaths_use_stochastic) deaths_multiplier else stop("stochastic scaling for deaths is not implemented yet.")
  
  df$yll <- df$yll * deaths_multiplier
  df$yld <- (df$dalys-df$yll) * cases_multiplier
  df$dalys <- df$yld + df$yll
  
  return(df)
}