## script to create the custom targeting table using as inputs a rds with DRC health zone locations
## and the DRC table csv with information on vaccination per health zone per year

## libraries

library(dplyr)
library(tidyr)
library(sf)
library(readr)
library(magrittr)
library(data.table)

## load data - your working directory should be 'gavi_vimc_cholera'

raw_table <- readr::read_csv("input_data/shapefiles/DRC_custom_shapefile/DRC_targeting_raw.csv") %>%
  dplyr::rename(NAME_2 = Health_zone) %>%
  dplyr::select(-gadm_match, -Notes)

locations <- readRDS("input_data/shapefiles/DRC_custom_shapefile/custom_shapefile.rds") 

## create targeting table

##left join locations with raw_table

raw_table_with_geometry <- dplyr::left_join(locations, raw_table, by = "NAME_2") %>%
  dplyr::mutate(actual_ocv1_fvp = 0, actual_ocv2_fvp = as.numeric(Vaccinated_people_2_rounds)/2) %>%
  dplyr::mutate(
   actual_prop_atleast_1dose_vaccinated = case_when(
    actual_ocv2_fvp > 0 ~ "0.965",
    is.na(actual_ocv2_fvp) ~ "0"
    )
   )

##add 0 vaccines when actual_ocv2_fvp is NA
raw_table_with_geometry[is.na(raw_table_with_geometry$actual_ocv2_fvp),"actual_ocv2_fvp"] <- 0

## final output
custom_targeting_table <- raw_table_with_geometry %>%
  dplyr::mutate(actual_ocv2_fvp = as.integer(actual_ocv2_fvp)) %>%
  dplyr::rename(vacc_year = Year) %>%
  dplyr::select(NAME_2, vacc_year, actual_prop_atleast_1dose_vaccinated, actual_ocv1_fvp, actual_ocv2_fvp)

##loop to add no vaccination years for each health zone
added_rows <- list()
vaccination_years <- c(2024, 2025, 2026) ## a list of possible vaccination years, each health zone gets vaccinated once
for (i in 1:nrow(custom_targeting_table)){
  ##get vaccination year
  year <- custom_targeting_table$vacc_year[i]
  ##if there is no vaccination for that health zone
  if (custom_targeting_table$actual_ocv2_fvp[i] == 0){
    custom_targeting_table$vacc_year[i] <- 2024 ## vacc_year will be NA so turn it to 2024
    other_year_1 <- custom_targeting_table[i,] %>% dplyr::mutate(vacc_year = 2025) ##the other two vaccination years
    other_year_2 <- custom_targeting_table[i,] %>% dplyr::mutate(vacc_year = 2026)
    added_rows[[i]] <- rbind(other_year_1, other_year_2)
  } else if (custom_targeting_table$actual_ocv2_fvp[i] > 0){
    ## for health zones with vaccination years
    no_vaccination_years <- vaccination_years[-which(vaccination_years == year)] ## get years without vaccination
    other_year_1 <- custom_targeting_table[i,] %>%
      dplyr::mutate(vacc_year = no_vaccination_years[1], actual_prop_atleast_1dose_vaccinated = 0,
                    actual_ocv1_fvp = 0, actual_ocv2_fvp = 0)
    other_year_2 <- custom_targeting_table[i,] %>%
      dplyr::mutate(vacc_year = no_vaccination_years[2], actual_prop_atleast_1dose_vaccinated = 0,
                    actual_ocv1_fvp = 0, actual_ocv2_fvp = 0)
    added_rows[[i]] <- rbind(other_year_1, other_year_2)
  }
}

## rbind added rows with no vaccination years to custom targeting table
added_rows <- sf::st_as_sf(data.table::rbindlist(added_rows)) ##add the no vaccination year rows to the targeting table
custom_targeting_table_final <- rbind(custom_targeting_table, added_rows) %>%
  dplyr::mutate(actual_prop_atleast_1dose_vaccinated = as.numeric(actual_prop_atleast_1dose_vaccinated))

##calculate target population for years with vaccination
custom_targeting_table_final <- custom_targeting_table_final %>%
  dplyr::mutate(pop_model = case_when(
    actual_ocv2_fvp == 0 ~ 0, ## for no vaccination years the target population is 0
    actual_ocv2_fvp != 0 ~ as.integer(actual_ocv2_fvp/actual_prop_atleast_1dose_vaccinated)
   )
  ) %>%
  dplyr::mutate(GID_0 = "COD") %>%
  dplyr::arrange(vacc_year)

## write custom targeting table to file

readr::write_csv(custom_targeting_table_final, file = "input_data/custom_targeting.csv")