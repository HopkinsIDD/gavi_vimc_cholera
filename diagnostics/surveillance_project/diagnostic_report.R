###################
# This R script will be turned into a .Rmd file after development
###################

# Parameters read into the script from the presets
params <- new.env()
params$surveillance_project_directory <- "/home/kaiyuezou/VIMC_Model/surveillance_project/gavi_vimc_cholera"
params$country <- c("COD") #country code specified or all the countries will be used for the diagnostic report
params$vac_incid_thresholds <- c(1/1000)
params$no_vaccination_surveillance_scenario <- "district-estimate" #specify which surveillance scenario is used for the no vaccination scenario

# Make a cache that saves all the tables and some important parameters 
cache <- new.env()

# Function that helps read in the config
get_filenames(cache = cache, surveillance_project_directory = params$surveillance_project_directory, 
    pre_country = params$country, pre_vac_incid_thresholds = params$vac_incid_thresholds, 
    no_vaccination_surveillance_scenario = params$no_vaccination_surveillance_scenario)

