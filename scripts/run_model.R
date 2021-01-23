# roxygen2::roxygenise("packages/ocvImpact")
# install.packages("packages/ocvImpact", type = "source", repos = NULL)

### Set Error Handling
if (Sys.getenv("INTERACTIVE_RUN", FALSE)) {
  options(warn = 1, error = recover)
} else {
  options(
    warn = 1,
    error = function(...) {
      quit(..., status = 2)
    }
  )
}

#### Libraries
package_list <- c(
  "data.table",
  "dplyr",
  "exactextractr",
  "fasterize",
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

### Run options
option_list <- list(
  optparse::make_option(
    c("-c", "--config"),
    action = "store",
    type = "character",
    help = "Model run configuration file"
  )
)
opt <- optparse::OptionParser(option_list = option_list) %>% optparse::parse_args()

### Read config file parameters
config <- yaml::read_yaml(opt$config)

runname <- config$runname
country <- config$country
scenario <- config$scenario
nsamples <- config$nsamples
targeting <- config$vacc$targeting_strategy
nskipyears <- config$vacc$num_skip_years

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
  clean = TRUE,
  targeting_strat = targeting,
  num_skip_years = nskipyears
  )

#### Create stochastic output file for country
message(paste("Calculating stochastic output:", runname, country, scenario, nsamples, targeting))
stoch_template <- ocvImpact::export_country_stoch_template(
  mpathname,
  country,
  scenario,
  ropathname
  )


message(paste("End script:", runname, country, scenario, nsamples, targeting))

