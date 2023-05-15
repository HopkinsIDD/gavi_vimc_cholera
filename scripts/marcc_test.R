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

### If run interactively
Sys.setenv(R_LIBRARY_DIRECTORY='/home/kzou7/rlibs/gcm/4.0.2/gcc/9.3.0/')
Sys.setenv(RUN_ON_MARCC=TRUE); 
Sys.setenv(CHOLERA_CONFIG='/home/kzou7/data-aazman1/kzou7/gavi-modeling/gavi_vimc_cholera/configs/202110gavi-3/campaign-default/global-estimate/0.001/BEN_campaign-default_global-estimate_200.yml')

### Libraries -- using the consistent one 7/2021
if (Sys.getenv("RUN_ON_MARCC", FALSE)) {
  r_lib <- Sys.getenv("R_LIBRARY_DIRECTORY", FALSE)
  library(codetools, lib=r_lib)
  library(remotes, lib=r_lib)
  library(sp, lib=r_lib)
  library(classInt, lib=r_lib)
  library(sf, lib=r_lib)
  library(rgdal, lib=r_lib)
  library(GADMTools, lib=r_lib)
  library(raster, lib=r_lib)
  library(terra, lib=r_lib)
  library(coda, lib=r_lib)
  library(ape, lib=r_lib)
  library(storr, lib=r_lib)
  library(abind, lib=r_lib)
  library(stars, lib=r_lib)

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

  library(desc, lib=r_lib)
  library(pkgload, lib=r_lib)
  roxygen2::roxygenise('packages/ocvImpact')
  install.packages('packages/ocvImpact', type='source', repos = NULL, lib=r_lib)

  library('montagu', character.only = T, lib = r_lib)
  library('ocvImpact', character.only = T, lib = r_lib)
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

## Screening
if(config$optimize$ir_pre_screening){
  message(paste(" -- Stage 1 Screening for", country))
  cache <- new.env()
  cache$rawoutpath <- ropathname
  cache$config <- config
  cache$ir_pre_screening_pass <- check_table_screening(spathname, country, scenario, as.numeric(config$vacc$vac_incid_threshold), tolower(config$vacc$vac_admin_level))
}

## Needed parameters 
datapath = dpathname
modelpath = mpathname
rawoutpath = ropathname
clean = cln
redraw = redrawIncid

ve_direct = generate_pct_protect_function()
indirect_mult = generate_indirect_incidence_mult()
secular_trend_mult = function(a,b,c,d){return(a*b*c*d)}

save_intermediate_raster <- as.logical(config$optimize$save_intermediate_raster)
save_final_output_raster <- as.logical(config$optimize$save_final_output_raster)
targeting_strategy <- config$vacc$targeting_strategy
vac_incid_threshold <- as.numeric(config$vacc$vac_incid_threshold)
vac_unconstrained <- as.logical(config$vacc$vac_unconstrained)
vac_admin_level <- tolower(config$vacc$vac_admin_level)
if(vac_admin_level == "both"){ rc_targeted <- c("rc1", "rc2")
  }else{ rc_targeted <- c("rc1", "rc2")[grepl(stringr::str_extract(vac_admin_level, "[0-9]{1}"), c("rc1", "rc2"))] }                       
vac_coverage <- as.numeric(config$vacc$vac_coverage)
surveillance_scenario <- config$surveillance_scenario$surveillance_scenario
vac_interval <- as.numeric(config$vacc$vac_interval)
sim_start_year <- as.numeric(config$vacc$sim_start_year)
vac_start_year <- as.numeric(config$vacc$vac_start_year)
vac_end_year <- as.numeric(config$vacc$vac_end_year)
sim_end_year <- as.numeric(config$vacc$sim_end_year)
num_skip_years <- as.numeric(config$vacc$num_skip_years)

incidence_rate_trend <- as.logical(config$setting$incidence_rate_trend)
use_country_incid_trend <- as.logical(config$incid$use_country_incid_trend)
outbreak_multiplier <- as.logical(config$setting$outbreak_multiplier)  
random_seed <- as.numeric(config$setting$random_seed)

## Testing functions 
shp0 <- load_shp0_by_country(datapath, country)
shp1 <- load_shp1_by_country(datapath, country)
shp2 <- load_shp2_by_country(datapath, country)
  
country_baseline <- ocvImpact::create_incid_raster(modelpath, datapath, country, nsamples, redraw, random_seed)
country_baseline_single_layer <- country_baseline[[1]]
pop_baseline <- ocvImpact::create_model_pop_raster(datapath, modelpath, country, 2022)
incid1 <- pop_weighted_admin_mean_incid(datapath, modelpath, incidence_rate_raster = country_baseline_single_layer, 
          pop_raster = pop_baseline, country, admin_shp = shp1)
incidence_rate_sum <- raster::extract(country_baseline_single_layer, shp1, method = 'bilinear', fun = sum, na.rm = TRUE)
incidence_rate_sum <- raster::extract(country_baseline, shp1, method = 'bilinear', fun = sum, na.rm = TRUE)

# incid_raster <- stars::st_as_stars(country_baseline)
# incid_vector <- nngeo::raster_extract(incid_raster, shp1, fun = sum, na.rm = TRUE)
