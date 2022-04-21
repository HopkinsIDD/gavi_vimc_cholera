surveillance_true_confirmation_rate <- function(datapath){
  omicron_dataset <- readr::read_csv(paste0(datapath, "/confirmation_rate/something.csv"))
  
  return(omicron_dataset$mean)
}
