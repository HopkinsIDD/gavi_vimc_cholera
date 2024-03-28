##this script creates a coverage table (proportion of target population vaccinated at admin 0 level) based on information inside 
##a custom-made targeting table (number of vaccinated people at admin 2 level)

## Load packages using Kaiyue's approach for other scripts in the repo
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
  "readxl",
  "tibble",
  "tidyr",
  "yaml",
  "ISOcodes"
)

for (package in package_list) {
  if (!require(package = package, character.only = T)) {
    install.packages(pkgs = package)
    library(package = package, character.only = T)
  }
  detach(pos = which(grepl(package, search())))
}

##load custom targeting table to calculate proportion vaccinated with one and two doses and get admin 0 population

targeting_table <- read_csv("input_data/custom_targeting.csv")  %>%
        group_by(vacc_year) %>%
        mutate(total_vaccinated_one_dose = sum(actual_ocv1_fvp)) %>% ## VP with one dose admin 0
        mutate(total_vaccinated_two_dose = sum(actual_ocv2_fvp)) %>% ## VP with two dose admin 0
        mutate(target = sum(pop_model)) %>% ## admin 0 population
        ungroup() %>%
        select(GID_0, vacc_year, total_vaccinated_one_dose, total_vaccinated_two_dose, target) %>%
        unique() %>%
        mutate(OCV1 = total_vaccinated_one_dose/target) %>% ## proportion vaccinated with one dose admin 0
        mutate(OCV2 = total_vaccinated_two_dose/target) %>% ## proportion vaccinated with two doses admin 0
        select(-total_vaccinated_one_dose, -total_vaccinated_two_dose) 

##format it and add columns to match the montagu coverage table format

coverage <- targeting_table %>%
  pivot_longer( cols = starts_with("OCV"), names_to = "vaccine", values_to = "coverage") %>%
  mutate(scenario = rep("cholera-ocv1-ocv2-default")) %>%
  mutate(gavi_support = rep("total")) %>%
  mutate(activity_type = rep("campaign")) %>%
  mutate(age_first = rep(1)) %>%
  mutate(age_last = rep(100)) %>%
  mutate(proportion_risk = rep('<NA>')) %>%  ##column is not used in cholera model, including it for consistency
  mutate(age_range_verbatim = rep('<NA>')) %>%
  mutate(gender = rep("both")) %>%
  mutate(set_name = ifelse(vaccine == 'OCV1', 'Cholera:OCV1,campaign,with:default', 'Cholera:OCV2,campaign,with:default')) %>%
  rename(country_code = GID_0) %>%
  rename(year = vacc_year) %>%
  ## get country name from its ISO code using the ISO_3166_1 from the ISOcodes package
  mutate(country = ISOcodes::ISO_3166_1$Name[which(ISO_3166_1$Alpha_3 == country_code)]) %>%
  ##reorder columns based on montagu coverage template
  select(scenario, set_name, vaccine, gavi_support, activity_type, country_code, country, year, age_first,
         age_last, age_range_verbatim, target, coverage, gender, proportion_risk) %>%
  arrange(vaccine) ##arrange rows by vaccine as in the montagu coverage template

##write to file
write_csv(coverage, "input_data/custom_coverage.csv")
