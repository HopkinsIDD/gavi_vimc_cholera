###################
# This R script will be turned into a .Rmd file after development
###################

# Parameters read into the script from the presets
params <- new.env()
params$surveillance_project_directory <- "/home/kaiyuezou/VIMC_Model/surveillance_project/gavi_vimc_cholera"
params$country <- c("COD") #country code specified or all the countries will be used for the diagnostic report
params$vac_incid_thresholds <- c(1/1000)
params$no_vaccination_surveillance_scenario <- c("district-estimate") #specify which surveillance scenario is used for the no vaccination scenario

# Packages 
library(dplyr)
library(tidyr)
library(ggplot2)

# Make a cache that saves all the tables and some important parameters 
cache <- new.env()

# Function that helps get the names of all files
get_filenames(cache = cache, surveillance_project_directory = params$surveillance_project_directory, 
    pre_country = params$country, pre_vac_incid_thresholds = params$vac_incid_thresholds, 
    no_vaccination_surveillance_scenario = params$no_vaccination_surveillance_scenario)

# Functions that aggregate the output across countries
combine_output(cache = cache, output_to_combine = c("target_table", "time_table"), save_final_output = TRUE)

# Function that returns cases plots (case type: true_case, observed_case, averted_tr_case, averted_ob_case) (cumulative_type: non_cumulative, district_level, country_level)
plt <- plot_cases(cache = cache, case_type = "true_case", threshold = 0.001, cumulative_type = "district_level")

# Function that returns efficacy plots (compare_type: side-by-side, subtraction) (cumulative_type: non_cumulative, across_year)
plt <- plot_efficacy(cache = cache, compare_type = "side-by-side", threshold = 0.001, cumulative_type = "non_cumulative")
