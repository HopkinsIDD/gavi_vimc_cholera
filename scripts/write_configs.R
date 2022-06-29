# roxygen2::roxygenise("packages/ocvImpact")
# install.packages("packages/ocvImpact", type = "source", repos = NULL)

#######Kaiyue Added on 7/21/2021#######
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
                  "yaml", 
                  "remotes", 
                  "rgeoboundaries"
                  )

for (package in package_list) {
  if (!require(package = package, character.only = T) & package != "rgeoboundaries") {
    install.packages(pkgs = package)
    library(package = package, character.only = T)
  }else if (package == "rgeoboundaries") {
    remotes::install_gitlab("dickoa/rgeoboundaries")
    if(!require(package = package, character.only = T)) {remotes::install_github("wmgeolab/rgeoboundaries")}
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



source("scripts/set_all_parameters.R")
#runname <- "201910gavi-5-config-test" #this is only for development use
cpathname <- file.path("configs", runname)
dir.create(cpathname, showWarnings = FALSE)

for(scn in scenarios){

  scnpathname <- file.path(cpathname, scn)
  dir.create(scnpathname, showWarnings = FALSE)

  for(surveillance_scenario in surveillance_scenarios){
    
    scnpathname <- file.path(cpathname, scn, surveillance_scenario)
    dir.create(scnpathname, showWarnings = FALSE)

    for(vac_incid_threshold in vac_incid_thresholds){
      
      scnpathname <- file.path(cpathname, scn, surveillance_scenario, vac_incid_threshold)
      dir.create(scnpathname, showWarnings = FALSE)

      pars <- tidyr::expand_grid(runname = runname, scenario = scn, country = countries, targeting = targeting_strategy, nsamples = num_samples, nskipyear = num_skip_years, clean = clean_outputs, redrawIncid = clean_incid) 
      pars$use_country_incid_trend <- use_country_incid_trend
      pars$incidence_rate_trend <- incidence_rate_trend
      pars$outbreak_multiplier <- outbreak_multiplier          
      pars$random_seed <- random_seed

      # the followings are specific to the surveillance project
      pars$save_intermediate_raster <- save_intermediate_raster
      pars$save_final_output_raster <- save_final_output_raster
      pars$vac_incid_threshold <- vac_incid_threshold
      pars$vac_unconstrained <- vac_unconstrained  
      pars$vac_admin_level <- vac_admin_level
      pars$vac_coverage <- vac_coverage
      pars$surveillance_scenario <- surveillance_scenario 
      pars$vac_interval <- vac_interval
      pars$sim_start_year <- sim_start_year
      pars$vac_start_year <- vac_start_year
      pars$vac_end_year <- vac_end_year
      pars$sim_end_year <- sim_end_year           

      lapply(1:nrow(pars), function(i){
        par_row <- pars[i,]
        ocvImpact::prepare_config(par_row, scnpathname)
      })
    }
  
  }  

}
