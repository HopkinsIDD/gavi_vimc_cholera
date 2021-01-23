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


#### Create paths
mpathname <- file.path("montagu", runname)
dpathname <- file.path("input_data")
ropathname <- file.path("output_raw", runname)
opathname <- file.path("output_final", runname)
dir.create(mpathname, showWarnings = FALSE)
dir.create(dpathname, showWarnings = FALSE)
dir.create(ropathname, showWarnings = FALSE)
dir.create(opathname, showWarnings = FALSE)

#### Crop incid raster by country
message(paste("Cropping incidence raster:", runname, country, scenario, nsamples))
incid <- ocvImpact::create_incid_raster(
  dpathname,
  country,
  ropathname,
  nsamples,
  clean = TRUE
  )

rm(incid)
gc()

message(paste("End incid_crop script:", runname, country, scenario, nsamples))

