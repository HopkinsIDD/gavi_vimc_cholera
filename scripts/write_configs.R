# roxygen2::roxygenise("packages/ocvImpact")
# install.packages("packages/ocvImpact", type = "source", repos = NULL)

#######Kaiyue Added on 7/21/2021#######
#======Use other packages needed======#
chooseCRANmirror(ind = 77) #specify the mirror so that the packages can be successfully installed in the non-interactive way

package_list <- c(
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

###########Comment completed###########



source("scripts/set_all_parameters.R")

cpathname <- file.path("configs", runname)
dir.create(cpathname, showWarnings = FALSE)

for(scn in scenarios){

  scnpathname <- file.path(cpathname, scn)
  dir.create(scnpathname, showWarnings = FALSE)
  pars <- tidyr::expand_grid(runname = runname, scenario = scn, country = countries, targeting = targeting_strategy, nsamples = num_samples, nskipyear = num_skip_years, clean = clean_outputs, redrawIncid = clean_incid)

  lapply(1:nrow(pars), function(i){
    par_row <- pars[i,]
    ocvImpact::prepare_config(par_row, scnpathname)
  })

}



