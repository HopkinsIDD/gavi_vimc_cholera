# roxygen2::roxygenise("packages/ocvImpact")
# install.packages("packages/ocvImpact", type = "source", repos = NULL)

source("scripts/set_all_parameters.R") ## read in parameters
opathname <- file.path("output_final", runname)
mpathname <- file.path("montagu", runname)

for (i in 1:length(scenarios)){

  scn <- scenarios[i]

  ## Stochastic outputs
  stoch_fns <- list.files(opathname, pattern = paste0(scn, "_stoch.csv"))
  stoch_out <- purrr::map_dfr(1:length(stoch_fns), function(j){
    return(readr::read_csv(file.path(opathname, stoch_fns[j])))
  })
  stoch_final_fn <- paste0(opathname, "/", ocvImpact::import_templateFilename_prefix("stochastic", mpathname), scn, ".csv")
  message(paste("Write stochastic final output:", stoch_final_fn))
  readr::write_csv(stoch_out, stoch_final_fn)

  ## Central outputs
  central_out <- ocvImpact::export_central_template(stoch_out)
  central_final_fn <- paste0(opathname, "/", ocvImpact::import_templateFilename_prefix("central", mpathname), scn, ".csv")
  message(paste("Write central final output:", stoch_final_fn))
  readr::write_csv(central_out, central_final_fn)

  ## Parameter outputs
  par_fns <- list.files(opathname, pattern = paste0(scn, "_pars.csv"))
  par_out <- purrr::map_dfr(1:length(par_fns), function(k){
    return(readr::read_csv(file.path(opathname, par_fns[k])))
  })
  par_final_fn <- paste0(opathname, "/", ocvImpact::import_templateFilename_prefix("parameter", mpathname), scenario, ".csv")
  message(paste("Write parameter final output:", par_final_fn))
  readr::write_csv(par_out, par_final_fn)
}
