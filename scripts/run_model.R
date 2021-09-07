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

#### Create paths
mpathname <- file.path("montagu", runname)
dpathname <- file.path("input_data")
ropathname <- file.path("output_raw", runname)
opathname <- file.path("output_final", runname)
dir.create(mpathname, showWarnings = FALSE)
dir.create(dpathname, showWarnings = FALSE)
dir.create(ropathname, showWarnings = FALSE)
dir.create(opathname, showWarnings = FALSE)

#### Run model
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

#### Create stochastic output file for country
message(paste("Calculating stochastic output:", runname, country, scenario, nsamples, targeting))
stoch_template <- ocvImpact::export_country_stoch_template(
  mpathname,
  country,
  scenario,
  ropathname,
  opathname
  )


message(paste("End script:", runname, country, scenario, nsamples, targeting))

