# This script is used to identify countries that are included in the ocv targeting plans and generate the corresponding targeting and coverage input files for vimc pipeline runs.
# ---- investment case study ---- #
# Helper functions ----
## identify_targeted_countries ----
identify_targeted_countries <- function(df, 
                                        coverage = 0.965, 
                                        rank_col = "rank_affected_pop",
                                        max_total_doses = 10 * 1e6,
                                        n_doses = 2) {
  required_cols <- c("country", "shapeName", "pop", rank_col)
  if (!all(required_cols %in% colnames(df))) {
    stop("Missing required columns in the input dataframe.")
  }
  
  df <- df %>%
    mutate(required_doses = pop * coverage * n_doses) %>%
    arrange(.data[[rank_col]]) %>%
    mutate(cumulative_doses = cumsum(required_doses))
  
  exceed_idx <- which(df$cumulative_doses > max_total_doses)[1]
  
  if (is.na(exceed_idx)) {
    df_targeting <- NULL
    message('No countries were identified')
  } else {
    df[exceed_idx, "required_doses"] <- max_total_doses - sum(df$required_doses[1:(exceed_idx - 1)])
    df$cumulative_doses <- cumsum(df$required_doses)
    
    df_targeting <- df[1:exceed_idx, ] %>%
      rename(GID_0 = country, GID_2 = shapeName)
  }
  
  return(df_targeting)
}

## create_custom_shapefiles ----
input_shapefiles <- readRDS("ocv_targeting_by_adm.rds")
create_custom_shapefiles <- function(df_targeted, input_shapefiles){
  library(tidyverse)
  library(sf)
  
  custom_shapefiles <- input_shapefiles %>% 
    filter(admin_level == "ADM2") %>% 
    filter(country %in% df_targeted$GID_0) %>% 
    rename(GID_0 = country, GID_2 = shapeName) %>% 
    select(GID_0,GID_2,admin_level,geom) %>% 
    group_by(GID_0,GID_2,admin_level,geom) %>% 
    slice(1) %>% 
    ungroup()
  
  return(custom_shapefiles)
}
custom_shapefiles <- create_custom_shapefiles(df_targeted,input_shapefiles)

## create_custom_country_shapefiles ----
create_custom_country_shapefiles <- function(df_targeted, input_shapefiles){
  library(tidyverse)
  library(sf)
  
  custom_country_shapefiles <- input_shapefiles %>% 
    filter(admin_level == "ADM0") %>% 
    filter(country %in% df_targeted$GID_0) %>% 
    rename(GID_0 = country, GID_2 = shapeName) %>% 
    select(GID_0,GID_2,admin_level,geom) %>% 
    group_by(GID_0,GID_2,admin_level,geom) %>% 
    slice(1) %>% 
    ungroup()

  return(custom_country_shapefiles)
}

custom_country_shapefiles <- create_custom_country_shapefiles(df_targeted,input_shapefiles)

## create_custom_targeting_file ----
create_custom_targeting_file <- function(df_targeted_raw, coverage = 0.965,custom_shapefiles,target_year = 2026, max_total_doses = 10 * 1e6,n_doses = 2){
  library(sf)
  library(tidyverse)
  
  custom_targeting_table  <- df_targeted_raw %>% 
    select(GID_0,GID_2,pop,cumulative_doses,required_doses) %>% 
    right_join(custom_shapefiles,by = c("GID_0","GID_2")) %>% 
    mutate(actual_ocv1_fvp = 0, 
                  actual_ocv2_fvp = as.numeric(required_doses/n_doses)) %>%
    mutate(
      actual_prop_atleast_1dose_vaccinated = case_when(
        actual_ocv2_fvp > 0 ~  actual_ocv2_fvp/pop,
        is.na(actual_ocv2_fvp) ~ 0
      )
    ) %>% 
    mutate(
      actual_ocv2_fvp = case_when(
        is.na(actual_ocv2_fvp) ~ 0,
        TRUE ~ actual_ocv2_fvp
      )
    ) %>% 
    mutate(vacc_year = target_year,
           actual_prop_atleast_1dose_vaccinated = as.numeric(actual_prop_atleast_1dose_vaccinated),
           pop_model = case_when(
             is.na(pop) ~ 0,
             TRUE ~ pop
           )) %>% 
    st_as_sf() %>% 
    select(GID_0,GID_2, vacc_year, pop_model,actual_prop_atleast_1dose_vaccinated, actual_ocv1_fvp, actual_ocv2_fvp)
  
  return(custom_targeting_table)

}

## create_custom_coverage_file ----
create_custom_coverage_file <- function(custom_targeting){
  library(tidyr)
  library(readr)
  library(ISOcodes) 
  library(dplyr)
  
  targeting_table <- custom_targeting %>%
    as.data.frame() %>%
    dplyr::group_by(vacc_year,GID_0) %>%
    dplyr::mutate(total_vaccinated_one_dose = sum(actual_ocv1_fvp)) %>% 
    dplyr::mutate(total_vaccinated_two_dose = sum(actual_ocv2_fvp)) %>% 
    dplyr::mutate(target = sum(pop_model))  %>% 
    dplyr::ungroup() %>%
    dplyr::select(GID_0, vacc_year, total_vaccinated_one_dose, total_vaccinated_two_dose, target) %>%
    unique() %>%
    dplyr::mutate(OCV1 = total_vaccinated_one_dose/target) %>% ## proportion vaccinated with one dose admin 0
    dplyr::mutate(OCV2 = total_vaccinated_two_dose/target) %>% ## proportion vaccinated with two doses admin 0
    dplyr::select(-total_vaccinated_one_dose, -total_vaccinated_two_dose) 
  
  custom_coverage <- targeting_table %>%
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
    rowwise() %>% 
    dplyr::mutate(country = ISOcodes::ISO_3166_1$Name[which(ISO_3166_1$Alpha_3 == country_code)]) %>%
    dplyr::select(scenario, set_name, vaccine, gavi_support, activity_type, country_code, country, year, age_first,
                  age_last, age_range_verbatim, target, coverage, gender, proportion_risk) %>%
    dplyr::arrange(vaccine)

  return(custom_coverage)
}

custom_coverage <- create_custom_coverage_file(custom_targeting)

## generate_investment_case_study_input_files ----
generate_investment_case_study_input_files <- function(df, 
                                                       max_total_doses, 
                                                       n_doses = 2,
                                                       coverage = 0.965,
                                                       rank_col = "rank_affected_pop",
                                                       input_shapefiles,
                                                       target_year = 2026,
                                                       recreate_file = TRUE) {
  library(dplyr)
  library(glue)
  library(stringr)
  
  df_targeted <- identify_targeted_countries(
    df = df,
    max_total_doses = max_total_doses,
    coverage = coverage,
    rank_col = rank_col
  )
  
  method_name <- str_remove(rank_col, "^rank_")
  dose_millions <- formatC(max_total_doses / 1e6, format = "f", digits = 0)
  
  default_modeled_countries_filename <- glue("ICS_{method_name}_{dose_millions}M_default_modeled_countries.csv")
  default_modeled_countries <- df_targeted %>% distinct(GID_0)
  if (recreate_file) {
    write.csv(default_modeled_countries, default_modeled_countries_filename, row.names = FALSE)
  } else {
    message(default_modeled_countries_filename, " has been found.")
  }
  
  custom_shapefiles <- create_custom_shapefiles(df_targeted, input_shapefiles)
  custom_country_shapefiles <- create_custom_country_shapefiles(df_targeted, input_shapefiles)
  custom_targeting <- create_custom_targeting_file(
    df_targeted_raw = df_targeted,
    custom_shapefiles = custom_shapefiles,
    target_year = target_year,
    coverage = coverage,
    max_total_doses = max_total_doses
  )
  custom_coverage <- create_custom_coverage_file(custom_targeting)
  
  if (recreate_file) {
    saveRDS(custom_shapefiles, glue("ICS_{method_name}_{dose_millions}M_custom_shapefiles.rds"))
    saveRDS(custom_country_shapefiles, glue("ICS_{method_name}_{dose_millions}M_country_shapefiles.rds"))
    saveRDS(custom_targeting, glue("ICS_{method_name}_{dose_millions}M_targeting_{target_year}.rds"))
    write.csv(custom_coverage, glue("ICS_{method_name}_{dose_millions}M_coverage_{target_year}.csv"), row.names = FALSE)
  }
  
  countries <- unique(df_targeted$GID_0)
  for (ctry in countries) {
    ctry_clean <- str_replace_all(ctry, "[^A-Za-z0-9]+", "_")
    
    ctry_shapefile <- custom_shapefiles %>% filter(GID_0 == ctry)
    ctry_country_shapefile <- custom_country_shapefiles %>% filter(GID_0 == ctry)
    ctry_targeting <- custom_targeting %>% filter(GID_0 == ctry)
    ctry_coverage <- custom_coverage %>% filter(country_code == ctry)
    
    shapefile_name <- glue("ICS_{ctry_clean}_{method_name}_{dose_millions}M_custom_shapefiles.rds")
    countryfile_name <- glue("ICS_{ctry_clean}_{method_name}_{dose_millions}M_country_shapefiles.rds")
    targeting_name <- glue("ICS_{ctry_clean}_{method_name}_{dose_millions}M_targeting_{target_year}.rds")
    coverage_name <- glue("ICS_{ctry_clean}_{method_name}_{dose_millions}M_coverage_{target_year}.csv")
    
    if (recreate_file) {
      saveRDS(ctry_shapefile, shapefile_name)
      saveRDS(ctry_country_shapefile, countryfile_name)
      saveRDS(ctry_targeting, targeting_name)
      write.csv(ctry_coverage, coverage_name, row.names = FALSE)
    } else {
      message("Files for ", ctry, " already exist.")
    }
  }
}


# Strategy 1: based on affected population  (2016_2020_MAI * annual_pop) ----
df <- read.csv("targeting_by_affected_pop.csv")
generate_investment_case_study_input_files(
  df,
  max_total_doses = 10 * 1e6, 
  coverage = 0.965,
  rank_col = "rank_affected_pop",
  input_shapefiles,
  target_year = 2026,
  recreate_file = T,
  n_doses = 2
)

# Strategy 2: based on prosective deaths  (MAI * new_who_cfr * pop) ----
df <- read.csv("targeting_by_prospect_deaths.csv")
generate_investment_case_study_input_files(
  df,
  max_total_doses = 10 * 1e6, 
  coverage = 0.965,
  rank_col = "rank_prospect_deaths",
  input_shapefiles =input_shapefiles,
  target_year = 2026,
  recreate_file = T
)

# Strategy 3: based on MAI ----
df <- read.csv("targeting_by_MAI.csv")
generate_investment_case_study_input_files(
  df,
  max_total_doses = 10 * 1e6, 
  coverage = 0.965,
  rank_col = "rank_MAI",
  input_shapefiles =input_shapefiles,
  target_year = 2026,
  recreate_file = T
)