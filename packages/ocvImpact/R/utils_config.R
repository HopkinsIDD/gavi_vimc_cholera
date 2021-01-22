
#' @title prepare_config
#' @name prepare_config
#' @description Automate generation of model configuration files from a dataframe of parameters. Generates one config per row in the p dataframe.
#'
#' @param p dataframe with config parameters, one config per row
#' @param configpath path to config files
#' @return Message returning name of yaml config file 
#' @export
prepare_config <- function(p, configpath){

  config_name <- paste0(configpath, "/", paste(p$country, p$scenario, p$nsamples, sep = "_"), ".yml")
  sink(file = config_name)

  cat(paste0(
    "runname: '", p$runname, "'\n",
    "country: '", p$country, "'\n",
    "scenario: '", p$scenario, "'\n",
    "num_samples: ", p$nsamples, "\n",
    "vacc:\n",
    "  targeting_strategy: ", p$targeting, "\n",
    "  num_skip_years: ", p$nskipyear, "\n"
  ))

  sink()
  return(message(paste(config_name, "written")))

}


