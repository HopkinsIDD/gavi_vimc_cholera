
### Set Error Handling
if (Sys.getenv("INTERACTIVE_RUN", FALSE)) {
  options(warn = 1, error = recover)
} else {
  options(
    warn = 1
  #  error = function(...) {
  #    quit(..., status = 2)
  #  }
  )
}

### Formal run checks (only for MARCC for now)
if (as.logical(Sys.getenv("RUN_ON_MARCC",FALSE))) {

  r_lib <- Sys.getenv("R_LIBRARY_DIRECTORY", FALSE)
  library(gert, lib=r_lib)

  #### Check local changes 
  if((nrow(gert::git_status(repo=getwd())) != 0)){
    mod_fns <- gert::git_status(repo=getwd())$file
    checklist <- c("^configs/", "^packages/ocvImpact/R/", "^scripts/")
    ignorelist <- c(".DS_Store$", "run_model.R$", "run_country_incid_crop.R$")
    for (itm in checklist) {
      for (ign in ignorelist){
        if(any(grepl(itm, mod_fns))){
          if(mod_fns[grepl(itm, mod_fns)] != "scripts/montagu_handle.R" & !all( grepl(ign, mod_fns[grepl(itm, mod_fns)]) )){
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



### Libraries -- depending on which server the model is being run on
if (Sys.getenv("RUN_ON_MARCC", FALSE)) {
  r_lib <- Sys.getenv("R_LIBRARY_DIRECTORY", FALSE)
  library(remotes, lib=r_lib)
  library(sp, lib=r_lib)
  library(classInt, lib=r_lib)
  library(sf, lib=r_lib)
  library(rgdal, lib=r_lib)
  library(GADMTools, lib=r_lib)
  library(raster, lib=r_lib)
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
  library(readxl)
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
#    "geodata", 
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
    "truncnorm", 
    "MCMCglmm"
  )

  for (package in package_list) {
    if (!require(package = package, character.only = T)) {
      install.packages(pkgs = package, lib.loc=Sys.getenv("R_LIBRARY_DIRECTORY", NULL))
      library(package = package, character.only = T, lib.loc=Sys.getenv("R_LIBRARY_DIRECTORY", NULL))
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
#    roxygen2::roxygenise("packages/ocvImpact")
#    install.packages("packages/ocvImpact", type = "source", repos = NULL)
    library('ocvImpact', character.only = T)
  }

  #======For the convenience of debugging======#
  ###These a few lines can be deleted safely after the model can run smoothly on the server. 
  library(raster)
  roxygen2::roxygenise("packages/ocvImpact")
  install.packages("packages/ocvImpact", type = "source", repos = NULL)
  library('ocvImpact', character.only = T)
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

#### Create paths
mpathname <- file.path("montagu", runname)
dpathname <- file.path("input_data")
opathname <- file.path("output_final", runname)
dir.create(mpathname, showWarnings = FALSE)
dir.create(dpathname, showWarnings = FALSE)
dir.create(opathname, showWarnings = FALSE)

#### Crop incid raster by country
#######Kaiyue Added on 7/14/2021####### --- add mpathname as input
message(paste("Cropping incidence raster:", runname, country, scenario, nsamples))
### Because the incidence raster data for IND and BGD are pre-made, skip this step for them
### This is added on Nov. 25th 

incid <- ocvImpact::create_incid_raster(
    mpathname, 
    dpathname,
    country,
    nsamples,
    clean <- TRUE, 
    redraw = as.logical(config$incid$redraw),
    random_seed = as.numeric(config$setting$random_seed) #updated from "=" to "<-"
    )



rm(incid)
gc()

message(paste("End incid_crop script:", runname, country, scenario, nsamples))

