# roxygen2::roxygenise("packages/ocvImpact")
# install.packages("packages/ocvImpact", type = "source", repos = NULL)

#######Kaiyue Added on 8/15/2021#######
#======Use other packages needed======#
chooseCRANmirror(ind = 77) #specify the mirror so that the packages can be successfully installed in the non-interactive way

package_list <- c(
  "GADMTools", 
  "rgdal", 
  "drat", 
  "roxygen2", 
  "data.table",
  "dplyr",
  "exactextractr",
  "fasterize",
  "optparse",
  "purrr",
  "raster",
  "readr",
  "sf",
  "stringr",
  "tibble",
  "tidyr",
  "yaml"
)

for (package in package_list) {
  if (!require(package = package, character.only = T)) {
    install.packages(pkgs = package)
    library(package = package, character.only = T)
  }
  detach(pos = which(grepl(package, search())))
}

#======Initialize Montagu package======#
if (!require('montagu', character.only = T)) {
  drat:::add("vimc")
  install.packages('montagu')
  library('montagu', character.only = T)
}
source("scripts/montagu_handle.R")

#======Use the ocvImpact package======#
if (!require('ocvImpact', character.only = T)) {
  roxygen2::roxygenise("packages/ocvImpact")
  install.packages("packages/ocvImpact", type = "source", repos = NULL)
  library('ocvImpact', character.only = T)
}

#======For the convenience of debugging======#
###These a few lines can be deleted safely after the model can run smoothly on the server. 
library(raster)
roxygen2::roxygenise("packages/ocvImpact")
install.packages("packages/ocvImpact", type = "source", repos = NULL)
library('ocvImpact', character.only = T)

###########Comment completed###########

library(magrittr)
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
  message(paste("Write central final output:", central_final_fn))
  readr::write_csv(central_out, central_final_fn)

  ## Parameter outputs
  par_fns <- list.files(opathname, pattern = paste0(scn, "_pars.csv"))
  par_out <- purrr::map_dfr(1:length(par_fns), function(k){
    return(readr::read_csv(file.path(opathname, par_fns[k])))
  }) %>%
    tidyr::pivot_wider(names_from = country, names_sep = ":", values_from = c(aoi, incid_rate, cfr))
  colnames <- stringr::str_split(names(par_out), ":", simplify=TRUE)
  names(par_out) <- stringr::str_remove(paste(colnames[,2], colnames[,1], sep = ":"), pattern = "^:")
  
  par_final_fn <- paste0(opathname, "/", ocvImpact::import_templateFilename_prefix("parameter", mpathname), scn, ".csv")
  message(paste("Write parameter final output:", par_final_fn))
  readr::write_csv(par_out, par_final_fn)
}
