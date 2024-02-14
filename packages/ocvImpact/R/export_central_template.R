#' @name export_central_template
#' @title export_central_template
#' @description Export a central burden template file (mean) using any stochastic template file output
#' @param stoch_output stochastic outputs, as from `[export_country_stoch_template()]`
#' @importFrom magrittr %>%
#' @return dataframe
#' @export
export_central_template <- function(stoch_output){
  
  if ('yll' %in% colnames(stoch_output)){
    central_output <- dplyr::group_by(stoch_output, disease, year, age, country, country_name) %>%
      dplyr::summarise(cohort_size = mean(cohort_size),
                       cases = mean(cases),
                       deaths = mean(deaths),
                       dalys = mean(dalys),
                       yll = mean(yll)) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(country, age, year) %>%
      dplyr::select(disease, year, age, country, country_name, cohort_size, cases, deaths, dalys, yll)
  } else {
    central_output <- dplyr::group_by(stoch_output, disease, year, age, country, country_name) %>%
      dplyr::summarise(cohort_size = mean(cohort_size),
                       cases = mean(cases),
                       deaths = mean(deaths),
                       dalys = mean(dalys)) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(country, age, year) %>%
      dplyr::select(disease, year, age, country, country_name, cohort_size, cases, deaths, dalys)
  }
  return(central_output)
}
