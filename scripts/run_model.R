### Set Error Handling
if (Sys.getenv("INTERACTIVE_RUN", FALSE)) {
  options(warn = 1, error = recover)
} else {
  options(
    warn = 1
    # error = function(...) {
    #   quit(..., status = 2)
    # }
  )
}

### Formal run checks (only for MARCC for now)
if (as.logical(Sys.getenv("RUN_ON_MARCC",FALSE))) {

  r_lib <- Sys.getenv("R_LIBRARY_DIRECTORY", FALSE)
  skip_checks <- as.logical(Sys.getenv("TESTING_RUN", FALSE))
  print(paste0("Whether to skip checks: ", skip_checks))
  library(gert, lib=r_lib)

  #### Check local changes 
  if((nrow(gert::git_status(repo=getwd())) != 0) & !skip_checks){
    mod_fns <- gert::git_status(repo=getwd())$file
    checklist <- c("^configs/", "^packages/ocvImpact/R/", "^scripts/")
    ignorelist <- c(".DS_Store$")
    for (itm in checklist) {
      for (ign in ignorelist){
        if(any(grepl(itm, mod_fns))){
          if(any(mod_fns[grepl(itm, mod_fns)] != "scripts/montagu_handle.R") & !all( grepl(ign, mod_fns[grepl(itm, mod_fns)]) )){
            mod_fns <- mod_fns[mod_fns != "scripts/montagu_handle.R"]
            stop(paste0("There are local changes that will affect formal run, please check: ", 
                        mod_fns[grepl(itm, mod_fns)][!grepl(ign, mod_fns[grepl(itm, mod_fns)])], "\n"))
          }
        }
      }
    }
  }

  #### Check package versions 
  if(packageVersion("sf", lib=r_lib) != "1.0.8" 
    |packageVersion("GADMTools", lib=r_lib) != "3.9.1"
    |packageVersion("exactextractr", lib=r_lib) != "0.9.0"
    |packageVersion("raster", lib=r_lib) != "3.4.13" 
    |packageVersion("Rcpp", lib=r_lib) != "1.0.9" 
    |packageVersion("terra", lib=r_lib) != "1.4.22"){ 
    stop("The important R packages do not have the correct versions, please check. ")
  }

}



### Libraries -- using the consistent one 7/2021
if (Sys.getenv("RUN_ON_MARCC", FALSE)) {
  r_lib <- Sys.getenv("R_LIBRARY_DIRECTORY", FALSE)
  library(vctrs, lib=r_lib)
  library(lifecycle, lib=r_lib)
  # library(codetools, lib=r_lib)
  # library(remotes, lib=r_lib)
  library(sp, lib=r_lib)
  library(classInt, lib=r_lib)
  library(sf, lib=r_lib)
  library(rgdal, lib=r_lib)
  library(tidyverse)
  library(dplyr, lib=r_lib)
  library(GADMTools, lib=r_lib)
  library(raster, lib=r_lib)
  library(terra, lib=r_lib)
  library(coda, lib=r_lib)
  library(ape, lib=r_lib)
  library(storr, lib=r_lib)

  library(drat, lib=r_lib)
  library(roxygen2, lib=r_lib)
  library(data.table, lib=r_lib)
  library(exactextractr, lib=r_lib)
  library(fasterize, lib=r_lib)
  library(truncnorm, lib=r_lib)
  library(MCMCglmm, lib=r_lib)

  # library(tibble)
  # library(withr)
  # library(processx)
  # library(optparse)
  # library(yaml)
  # library(purrr)
  library(dplyr)
  # library(tidyverse)
  # library(stringi)
  # library(readr)
  # library(stringr)
  # library(tidyr)

  ###comment out the following lines when launch formal runs (no need to comment out if on MARCC)
  library(desc, lib=r_lib)
  library(pkgload, lib=r_lib)
  roxygen2::roxygenise('packages/ocvImpact')
  install.packages('packages/ocvImpact', type='source', repos = NULL, lib=r_lib)

  library(montagu, lib = r_lib)
  library(ocvImpact, lib = r_lib)
  source("scripts/montagu_handle.R")


} else {
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
    "truncnorm", 
    "MCMCglmm"
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
  ### The two lines below should be commented out when running the formal model (as it may cause issues when running multiple countries at the same time)
  roxygen2::roxygenise("packages/ocvImpact")
  install.packages("packages/ocvImpact", type = "source", repos = NULL)

  ##======Load certain packages that are used a lot======#
  library('ocvImpact', character.only = T)
  library(raster)
  library(dplyr)
}





#======================================== the start of the run ========================================#
### Run options
option_list <- list(
  optparse::make_option(
    c("-c", "--config"),
    action = "store",
    default = Sys.getenv("CHOLERA_CONFIG", "config.yml"),
    type = "character",
    help = "Model run configuration file"
  )
)
opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))



### Read config file parameters
config <- yaml::read_yaml(opt$config)

runname <- config$runname
country <- config$country
scenario <- config$scenario
nsamples <- config$incid$num_samples
redrawIncid <- config$incid$redraw
targeting <- config$vacc$targeting_strategy
nskipyears <- config$vacc$num_skip_years
cln <- config$clean
##number of doses, a parameter that only applies to the 202310gavi-4 touchstone
if(runname == "202310gavi-4"){
  ndoses <- config$ndoses
}

# random_seed <- as.numeric(config$setting$random_seed) #this one will be used through redrawing incidence rate raster and generating other rasters
# set.seed(random_seed)
# rm(random_seed)

#### Create paths
mpathname <- file.path("montagu", runname)
dpathname <- file.path("input_data")
spathname <- file.path("input_data", "incidence")
ropathname <- file.path("output_raw", runname)
opathname <- file.path("output_final", runname)
dir.create(mpathname, showWarnings = FALSE)
dir.create(dpathname, showWarnings = FALSE)
dir.create(spathname, showWarnings = FALSE)
dir.create(ropathname, showWarnings = FALSE)
##calam added to create separate paths for the one-dose and two-dose scenarios for the 202310gavi-4 touchstone
if (runname == "202310gavi-4"){
  dir.create(file.path(ropathname, scenario, ndoses), showWarnings = FALSE)
} else {
  dir.create(file.path(ropathname, scenario), showWarnings = FALSE)
}
##the addition ends here
dir.create(opathname, showWarnings = FALSE)

#### Run model -- where different projects diverge 
if(config$vacc$targeting_strategy == 'threshold_unconstrained'){
  ### The surveillance project
  cache <- new.env()
  cache$rawoutpath <- ropathname
  cache$config <- config

  ## Screening
  if(config$optimize$ir_pre_screening){
    message(paste(" -- Stage 1 Screening for", country))
    cache$ir_pre_screening_pass <- check_table_screening(spathname, country, scenario, as.numeric(config$vacc$vac_incid_threshold), tolower(config$vacc$vac_admin_level))
  }else{
    cache$ir_pre_screening_pass <- TRUE
    cache$novacc_campde_transfer <- FALSE
  }

  ## If the stage 1 screening process has been passed 
  message(paste(" --- Running Surveillance Project:", country, "for scenario combination:", 
                config$surveillance_scenario$surveillance_scenario, 
                scenario, as.numeric(config$vacc$vac_incid_threshold), 
                tolower(config$vacc$vac_admin_level), nsamples)) #needs to add more
  run_surveillance_scenario( 
    datapath = dpathname,
    modelpath = mpathname,
    country,
    scenario,
    rawoutpath = ropathname,
    nsamples,
    clean = cln,
    redraw = redrawIncid,
    config = config #all new parameters that would be used will be read in as one config file 
  ) 


}else if(config$vacc$targeting_strategy != 'threshold_unconstrained' & runname == "202310gavi-4"){
  ### The VIMC project for the 202310gavi-4 touchstone with one dose and two dose scenarios
  message(paste("Running:", runname, country, scenario, nsamples, targeting, ndoses))
  expected_cases <- ocvImpact::run_country_scenario(
    dpathname,
    mpathname,
    country,
    scenario,
    ropathname,
    nsamples,
    num_doses = ndoses,
    clean = cln,
    redraw = redrawIncid,
    targeting_strat = targeting,
    num_skip_years = nskipyears
  )
  
  ## Create stochastic output file for country
  message(paste("Calculating stochastic output:", runname, country, scenario, nsamples, targeting, ndoses))
  stoch_template <- ocvImpact::export_country_stoch_template(
    mpathname,
    country,
    scenario,
    ropathname,
    opathname
  )
  
  
  message(paste("End script:", runname, country, scenario, nsamples, targeting, ndoses))
  }else{
  ### The VIMC project 
  message(paste("Running:", runname, country, scenario, nsamples, targeting))
  expected_cases <- ocvImpact::run_country_scenario(
    dpathname,
    mpathname,
    country,
    scenario,
    ropathname,
    nsamples,
    clean = cln,
    redraw = redrawIncid,
    targeting_strat = targeting,
    num_skip_years = nskipyears
    )

  ## Create stochastic output file for country
  message(paste("Calculating stochastic output:", runname, country, scenario, nsamples, targeting))
  stoch_template <- ocvImpact::export_country_stoch_template(
    mpathname,
    country,
    scenario,
    ropathname,
    opathname
    )


  message(paste("End script:", runname, country, scenario, nsamples, targeting))
}

# rm(list = ls()) #temp *******************************************************************
# gc() #temp ******************************************************************************
