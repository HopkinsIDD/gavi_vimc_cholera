# roxygen2::roxygenise("packages/ocvImpact")
# install.packages("packages/ocvImpact", type = "source", repos = NULL)

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

#### Libraries -- using the consistent one 7/2021
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

###########Comment completed###########



### Run options
option_list <- list(
  optparse::make_option(
    c("-c", "--config"),
    action = "store",
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
dir.create(file.path(ropathname, scenario), showWarnings = FALSE)
dir.create(opathname, showWarnings = FALSE)

#### Run model -- where different projects diverge 
if(config$vacc$targeting_strategy == 'threshold_unconstrained'){
  ### The surveillance project
  
  ## Screening
  if(config$optimize$ir_pre_screening){
    message(paste(" -- Stage 1 Screening for", country))
    cache <- new.env()
    cache$rawoutpath <- ropathname
    cache$config <- config
    cache$ir_pre_screening_pass <- check_table_screening(spathname, country, scenario, as.numeric(config$vacc$vac_incid_threshold), tolower(config$vacc$vac_admin_level))
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
