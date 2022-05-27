#' @name surveillance_true_confirmation_rate
#' @title surveillance_true_confirmation_rate
#' @description surveillance_true_confirmation_rate
#' @param datapath path to input data
#' @return 
#' @export
#' @include
surveillance_true_confirmation_rate <- function(datapath){
  
  omicron_dataset <- readr::read_csv(paste0(datapath, "/confirmation_rate/parameters.csv"))
  
  return(omicron_dataset$mean)
}
