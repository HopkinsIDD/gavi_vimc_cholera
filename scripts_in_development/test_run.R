## This script tests the added functions with COD
library(sf)
library(dplyr)
setwd("gavi_vimc_cholera")

# load existing package
if (!require('ocvImpact', character.only = T)) {
   roxygen2::roxygenise("packages/ocvImpact")
   install.packages("packages/ocvImpact", type = "source", repos = NULL)
   library('ocvImpact', character.only = T)}

# set the parameters and montagu handle
source("scripts/set_all_parameters.R")
source("scripts/montagu_handle.R")
datapath <- "input_data"
modelpath <- "montagu/202110gavi-3"
rawoutpath <- "output_raw/202110gavi-3"
country <- "COD"
scenario = "campaign-default"

# source newly added functions (might need to change path)
source("~/vimc/draft/code/update_rc.R")
source("~/vimc/draft/code/update_input_rasterStack.R")
source("~/vimc/draft/code/update_sus_rasterStack.R")
source("~/vimc/draft/code/update_ec.R")
source("~/vimc/draft/code/update_incidence.R")
source("~/vimc/draft/code/run_one_country.R")



# test with COD (campaigns from 2018 to 2022)
run_one_country(datapath = datapath, country = "COD", campaign_cov = 0.8, modelpath, rawoutpath,
                baseline_year = 2018, end_year = 2022, nsamples = 30, redraw = FALSE)
