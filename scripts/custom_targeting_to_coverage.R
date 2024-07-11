##this script creates a coverage table (proportion of target population vaccinated at admin 0 level) based on information inside 
##a custom-made targeting table (number of vaccinated people at admin 2 level)

## Load libraries
chooseCRANmirror(ind = 77) #specify the mirror so that the packages can be successfully installed in the non-interactive way

library(tidyr)
library(readr)
library(ISOcodes) 
library(dplyr)

## load custom targeting table to calculate proportion vaccinated with one and two doses and get admin 0 population

## the custom targeting table should have one row per admin 2 unit and include the following columns:
## vacc_year: numeric, year of vaccination
## actual_ocv1_fvp: numeric, number of people who received one dose of OCV
## actual_ocv2_fvp: numeric, number of people who received two doses of OCV
## pop_model: numeric, population size per admin 2 unit
## GID_0: character, ISO 3166-1 alpha-3 code for the country

## modify the filepath below to specify the name of the custom targeting table and the directory it's in

targeting_table <- readRDS("input_data/drc_custom_targeting_2024_2026.rds")  %>%
  as.data.frame() %>%
  dplyr::group_by(vacc_year) %>%
  dplyr::mutate(total_vaccinated_one_dose = sum(actual_ocv1_fvp)) %>% ## VP with one dose admin 0
  dplyr::mutate(total_vaccinated_two_dose = sum(actual_ocv2_fvp)) %>% ## VP with two dose admin 0
  dplyr::mutate(target = sum(pop_model)) %>% ## admin 0 population
  dplyr::ungroup() %>%
  dplyr::select(GID_0, vacc_year, total_vaccinated_one_dose, total_vaccinated_two_dose, target) %>%
  unique() %>%
  dplyr::mutate(OCV1 = total_vaccinated_one_dose/target) %>% ## proportion vaccinated with one dose admin 0
  dplyr::mutate(OCV2 = total_vaccinated_two_dose/target) %>% ## proportion vaccinated with two doses admin 0
  dplyr::select(-total_vaccinated_one_dose, -total_vaccinated_two_dose) 

##format it and add columns to match the montagu coverage table format

coverage <- targeting_table %>%
  tidyr::pivot_longer( cols = starts_with("OCV"), names_to = "vaccine", values_to = "coverage") %>%
  dplyr::mutate(scenario = rep("cholera-ocv1-ocv2-default")) %>%
  dplyr::mutate(gavi_support = rep("total")) %>%
  dplyr::mutate(activity_type = rep("campaign")) %>%
  dplyr::mutate(age_first = rep(1)) %>%
  dplyr::mutate(age_last = rep(100)) %>%
  dplyr::mutate(proportion_risk = rep('<NA>')) %>%  ##column is not used in cholera model, including it for consistency
  dplyr::mutate(age_range_verbatim = rep('<NA>')) %>%
  dplyr::mutate(gender = rep("both")) %>%
  dplyr::mutate(set_name = ifelse(vaccine == 'OCV1', 'Cholera:OCV1,campaign,with:default', 'Cholera:OCV2,campaign,with:default')) %>%
  dplyr::rename(country_code = GID_0) %>%
  dplyr::rename(year = vacc_year) %>%
  ## get country name from its ISO code using the ISO_3166_1 from the ISOcodes package
  dplyr::mutate(country = ISOcodes::ISO_3166_1$Name[which(ISO_3166_1$Alpha_3 == country_code)]) %>%
  ##reorder columns based on montagu coverage template
  dplyr::select(scenario, set_name, vaccine, gavi_support, activity_type, country_code, country, year, age_first,
                age_last, age_range_verbatim, target, coverage, gender, proportion_risk) %>%
  dplyr::arrange(vaccine) ##arrange rows by vaccine as in the montagu coverage template

##write to file
readr::write_csv(coverage, "input_data/drc_custom_coverage_2024_2026.csv")
