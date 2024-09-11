## code to prepare `DATASET` dataset goes here

## This script creates package data from existing files inside the input_data directory

#### load the necessary files from the input_data directory as objects ####

## Parameters for beta distribution that fits the confirmation rate distribution, Used for the surveillance project
surveillance_confirmation_rate_parameters <- readr::read_csv("../../../input_data/confirmation_rate/parameters.csv")

## Country name, country code, and incidence rate for 47 VIMC countries for which we do not have modeled raster estimates
## These country-level estimates are projected to a grid and used to run the VIMC model for these countries. Last updated in 2020
default_country_incidence <- readr::read_csv("../../../input_data/incidence/VIMC-47-countries-for-cholera-modelling.csv")

##Reported data on country level suspected cases and deaths, used to calibrate incidence rate trends. source: WHO Annual Cholera Reports, last updated in 2020
who_annual_reports <- readr::read_csv("../../../input_data/incidence/who_case_repo_source.csv")

## Vaccine efficacy and effectiveness data from Bi et al. (2017) that was used to model VE decay in the 2019 and 2021 touchstones
## These data are longer used but this object was retained for backwards compatibility
ve_raw_2017 <- readr::read_csv("../../../input_data/ocv_ve_overtime.csv")

## Annual case fatality ratios by country. Calculated from case and death data in WHO Annual Cholera Reports. Last updated in 2019
who_cfrs <- readr::read_csv("../../../input_data/who_cfrs.csv")

## Modeled one dose vaccine effectiveness decay function fit from only observational studies, used in the 2023 touchstone, source: Xu and Azman 2024
ve_log_1dose_obs <- readRDS("../../../input_data/log1d_obs.rds")

## Modeled two dose vaccine effectiveness decay function fit from only observational studies, used in the 2023 touchstone, source: Xu and Azman 2024
ve_log_2dose_obs <- readRDS("../../../input_data/log2d_obs.rds")

## Disability weights used to calculate DALYs and YLLs in 2019, 2021, and 2023 touchstones, source: IHME Global Burden of Disease, 2019
disability_weights <- readxl::read_xlsx("../../../input_data/IHME_GBD_2019_DISABILITY_WEIGHTS_Y2020M010D15.XLSX")

## List of default countries that could be modeled in the VIMC core model, as called in scripts/set_all_parameters.R
## The default value is currently overwritten by another list of countries hard-coded in set_all_parameters 
default_countries <- readr::read_csv("../../../input_data/default_modeled_countries.csv")

## Custom coverage table used for the 2024 DRC vaccine impact case study 
## This represents an example file for scenarios with "custom" coverage settings (when custom$use_montagu_coverage: FALSE, custom$coverage_filename)
drc_custom_coverage_2024_2026 <- readr::read_csv("../../../input_data/drc_custom_coverage_2024_2026.csv")

## Custom targeting table used for the 2024 DRC vaccine impact case study
## This represents an example file for scenarios with a "custom" targeting strategy (custom$targeting_filename)
drc_custom_targeting_2024_2026 <- readRDS("../../../input_data/drc_custom_targeting_2024_2026.rds")

## Custom health zone shapefile used for the DRC vaccine impact case study
## This represents an example custom shapefile set (when custom$use_custom_shapefile: TRUE, custom$shapefile_filename)
drc_custom_healthzone_shapefile <- readRDS("../../../input_data/shapefiles/DRC_custom_shapefile/custom_shapefile.rds")

## Custom country shapefile used for the DRC vaccine impact case study
## This represents an example country-level custom shapefile set (when custom$use_custom_shapefile: TRUE, custom$country_shapefile_filename)
drc_custom_country_shapefile <- readRDS("../../../input_data/shapefiles/DRC_custom_shapefile/country_shapefile.rds")


#### Add them as package data using the "usethis" package ####

usethis::use_data(surveillance_confirmation_rate_parameters)
usethis::use_data(default_country_incidence)
usethis::use_data(who_annual_reports)
usethis::use_data(ve_raw_2017)
usethis::use_data(who_cfrs)
usethis::use_data(ve_log_1dose_obs)
usethis::use_data(ve_log_2dose_obs)
usethis::use_data(disability_weights)
usethis::use_data(default_countries)
usethis::use_data(drc_custom_coverage_2024_2026)
usethis::use_data(drc_custom_targeting_2024_2026)
usethis::use_data(drc_custom_healthzone_shapefile)
usethis::use_data(drc_custom_country_shapefile)


