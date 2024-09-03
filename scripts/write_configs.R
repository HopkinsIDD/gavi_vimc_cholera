
#######Kaiyue Added on 7/21/2021#######
#======Use other packages needed======#
chooseCRANmirror(ind = 77) #specify the mirror so that the packages can be successfully installed in the non-interactive way

package_list <- c(
                  "geodata", 
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
                  "readxl",
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


library('ocvImpact', character.only = T)

###########Comment completed###########



source("scripts/set_all_parameters.R")
#runname <- "201910gavi-5-config-test" #this is only for development use
cpathname <- file.path("configs", runname)
#dir.create(cpathname, showWarnings = FALSE)

for(scn in scenarios){

  scnpathname <- file.path(cpathname, scn)
  dir.create(scnpathname, showWarnings = FALSE)
  
  #for (dose in ndoses){

    for(surveillance_scenario in surveillance_scenarios){
    
    ##scnpathname <- file.path(cpathname, scn, surveillance_scenario)
    ##dir.create(scnpathname, showWarnings = FALSE)

      for(vac_incid_threshold in vac_incid_thresholds){
      

        # parameters that apply to both projects 
        pars <- tidyr::expand_grid(runname = runname, scenario = scn, country = countries, targeting = targeting_strategy, nsamples = num_samples, nskipyear = num_skip_years, clean = clean_outputs, redrawIncid = clean_incid) 
        pars$use_country_incid_trend <- use_country_incid_trend
        pars$incidence_rate_trend <- incidence_rate_trend
        pars$outbreak_multiplier <- outbreak_multiplier          
        pars$random_seed <- random_seed
        pars$campaign_cov <- campaign_cov
        pars$use_montagu_coverage <- use_montagu_coverage
        pars$use_custom_shapefile <- use_custom_shapefile
        
        if(targeting_strategy != "threshold_unconstrained" & runname != "202310gavi-4"){
          scnpathname <- file.path(cpathname, scn)
        }else if(targeting_strategy != "threshold_unconstrained" & runname == "202310gavi-4"){
            scnpathname <- file.path(cpathname, scn, incidence_rate_trend)
        }else if(targeting_strategy == "threshold_unconstrained"){
          scnpathname <- file.path(cpathname, scn, surveillance_scenario, vac_incid_threshold)
        }
        dir.create(scnpathname, showWarnings = FALSE)
        
        #parameters that apply only to the 202310gavi-4 touchstone with the one-dose and two-dose vaccination scenarios
        #pars$ndoses <- dose
        if (scn == "ocv1-default"){
          pars$ndoses <- "one"
        } else if (scn == "ocv1-ocv2-default"){
          pars$ndoses <- "two"
        } else {
          pars$ndoses <- "zero"
        }
        
        #parameters needed for the DRC Case study (202310gavi-4 touchstone)
        if(use_montagu_coverage == FALSE){
          pars$output_years <- output_years
          pars$custom_targeting_filename <- custom_targeting_filename
          pars$custom_coverage_filename <- custom_coverage_filename
        }
        
        if(use_custom_shapefile == TRUE){
          pars$custom_shapefile_filename <- custom_shapefile_filename
          pars$custom_country_shapefile_filename <- custom_country_shapefile_filename
        }
        
        # the followings are specific to the surveillance project
        pars$save_intermediate_raster <- save_intermediate_raster
        pars$save_final_output_raster <- save_final_output_raster
        pars$ir_pre_screening <- ir_pre_screening
        pars$vac_incid_threshold <- vac_incid_threshold
        pars$vac_unconstrained <- vac_unconstrained  
        pars$vac_admin_level <- vac_admin_level
        pars$vac_coverage <- vac_coverage
        pars$surveillance_scenario <- surveillance_scenario 
        pars$testing_sensitivity <- testing_sensitivity[[surveillance_scenario]]
        pars$vac_interval <- vac_interval
        pars$sim_start_year <- sim_start_year
        pars$vac_start_year <- vac_start_year
        pars$vac_end_year <- vac_end_year
        pars$sim_end_year <- sim_end_year   
        pars$use_mean_ir <- use_mean_ir
        pars$mean_ir_span <- mean_ir_span       

        lapply(1:nrow(pars), function(i){
          par_row <- pars[i,]
          ocvImpact::prepare_config(par_row, scnpathname)
        })
      }
  
    }  
 # }
}
