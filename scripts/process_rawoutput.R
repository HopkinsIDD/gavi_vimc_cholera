# Script to summarize the raw output to get intermediate output datasets
# the intermediate datasets are then to be used in the rmd to make the manuscript figures/tables

# set up =======================================================
library(dplyr, lib = r_lib)
library(stringr)
library(tidyr)
library(tidyverse)
library(sf, lib=r_lib)
library(sp, lib=r_lib)
library(classInt, lib=r_lib)
library(rgdal, lib=r_lib)
library(GADMTools, lib = r_lib)
library(raster, lib=r_lib)
library(exactextractr, lib=r_lib)
library(s2, lib=r_lib)
source("packages/ocvImpact/R/utils_new_op.R")
output_final_path <- "output_final/202302_survms"
all_countries <- c("AGO", "BDI", "BEN", "BFA", "CAF", "CIV", "CMR", "COD", "COG", "ETH", "GHA", "GIN", "GNB", "KEN", "LBR", "MDG", "MLI", "MOZ", "MRT", "MWI", "NAM", "NER", "NGA", "RWA", "SEN", "SLE", "SOM", "SSD", "TCD", "TGO", "TZA", "UGA", "ZAF", "ZMB", "ZWE")
# r_lib <- "/home/hxu70/rlibs/gcm/4.0.2/gcc/9.3.0/"

# check directory status =======================================================
if(!file.exists(paste0(output_final_path, "/", "target_table"))){
    dir.create(paste0(output_final_path, "/", "target_table"))
}else{
    if(length(dir(all.files=TRUE)) != 0){
        message("Directory output_final/202302_survms/target_table is not empty and target tables have already been made.")
    }
}


# OUTPUT 1: target_table =======================================================
# make target_table for each country and store them separately (one country per csv) into "output_final/202302_survms/target_table", this table combines all scenarios together for each country
for(country in all_countries){
    
    fn <- paste0("target_table_", country, ".csv")
    target_table  <- make_target_table(countries = country, thresholds = c("1e-04", "2e-04", "0.001"), output_raw_path = "output_raw/202302_survms")
    write.csv(target_table, paste0(output_final_path, "/target_table/", fn), row.names = F)
    
}


# OUTPUT 2: tp_eff_ac_byISO_medianCI table =======================================================
# make table of median and CI of targeted population, OCV efficiency and averted true cases by country
for(country in all_countries){
    
    # read in target table 
    df_tt <- read.csv(paste0(output_final_path, "/target_table/target_table_", country, ".csv" ))
    df_tt <- df_tt %>% mutate(true_averted_cases = no_vaccination_true_case - campaign_default_true_case)

    # added on 11/10/2022: remove targets smaller than 5*5 grid (averted case == 0 & is_target == 1)
    df_tt <- df_tt %>% filter(true_averted_cases != 0 | is_target != 1)
    
    # get targeted pop and efficiency table for this one country
    df_tp_eff <- combine_eff_table_median(countries = country, df_target = df_tt, admin_level = "both")
    
    # rbind all countries
    if(which(all_countries == country) == 1){df_tp_eff_byISO_medianCI <- df_tp_eff}
    if(which(all_countries == country) != 1){df_tp_eff_byISO_medianCI <- rbind(df_tp_eff_byISO_medianCI, df_tp_eff)}

}

if(!file.exists(paste0(output_final_path, "/intermediate_table"))){
    dir.create(paste0(output_final_path, "/intermediate_table"))
    write.csv(df_tp_eff_byISO_medianCI, paste0(output_final_path, "/intermediate_table/tp_eff_ac_byISO_medianCI.csv"), row.names = F)
}else{
    write.csv(df_tp_eff_byISO_medianCI, paste0(output_final_path, "/intermediate_table/tp_eff_ac_byISO_medianCI.csv"), row.names = F)
}


# OUTPUT 3: tp_ac_allISOs =======================================================
# make table of targeted population, efficiency and averted true cases for all countries combined together (summarize across Africa)
for(country in all_countries){ # stack target tables of all countries

    df_tt <- read.csv(paste0(output_final_path, "/target_table/target_table_", country, ".csv" ))
    if(which(all_countries == country) == 1){df_tt_allISOs <- df_tt}
    if(which(all_countries == country) != 1){df_tt_allISOs <- rbind(df_tt_allISOs, df_tt)}

}

df_tt_allISOs <- df_tt_allISOs %>%  mutate(true_averted_cases = no_vaccination_true_case - campaign_default_true_case)

# added on 11/10/2022: remove targets smaller than 5*5 grid (averted case == 0 & is_target == 1)
df_tt_alllISOs <- df_tt_allISOs %>% filter(true_averted_cases != 0 | is_target != 1)

df_tp <- df_tt_allISOs %>%
    #filter(is_target == 1) %>% # added on 10/19/2022 keeping runs even there is no targeting going on
    group_by(year, run_id, threshold, confirmation_lens, admin_level) %>%
    summarize(target_pop = sum(actual_fvp)) %>%
    group_by(run_id, threshold, confirmation_lens, admin_level) %>%
    mutate(target_pop_cumu = cumsum(target_pop)) %>%
    dplyr::select(-target_pop)

df_ac <- df_tt_allISOs %>%
    group_by(year, run_id, threshold, confirmation_lens, admin_level) %>%
    summarize(true_ac = sum(true_averted_cases)) %>%
    group_by(run_id, threshold, confirmation_lens, admin_level) %>%
    mutate(true_ac_cumu = cumsum(true_ac)) %>%
    dplyr::select(-true_ac)
    
tp_ac_allISOs <- df_tp %>% 
    full_join(df_ac, by = c("year", "run_id", "threshold", "confirmation_lens", "admin_level")) 

write.csv(tp_ac_allISOs, paste0(output_final_path, "/intermediate_table/tp_ac_allISOs.csv"), row.names = F)


# OUTPUT 4: tp_eff_ac_allISOs_medianCI =======================================================
# use tp_ac_allISOs (OUTPUT 3) to calculate median and CI of targeted population, efficiency, and averted true cases across all modeled countries
if(!file.exists(paste0(output_final_path, "/intermediate_table/tp_ac_allISOs.csv"))){
    message("Data file tp_ac_allISOs.csv is not available in the intermediate_table folder.")
}else{
    tp_ac_allISOs <- read.csv(paste0(output_final_path, "/intermediate_table/tp_ac_allISOs.csv"))
}

df_tp <- tp_ac_allISOs %>%
    group_by(year, threshold, confirmation_lens, admin_level) %>%
    summarize(tp_cumu_median = median(target_pop_cumu, na.rm = T),
              tp_cumu_lb = quantile(target_pop_cumu, 0.025, na.rm = T),
              tp_cumu_ub = quantile(target_pop_cumu, 0.975, na.rm = T)) %>%
    dplyr::select(year, threshold, confirmation_lens, admin_level, tp_cumu_lb, tp_cumu_median, tp_cumu_ub) %>%
    group_by(threshold, confirmation_lens, admin_level) %>%
    summarize(tp_cumu_lb = max(tp_cumu_lb, na.rm = T),
              tp_cumu_median = max(tp_cumu_median, na.rm = T),
              tp_cumu_ub = max(tp_cumu_ub, na.rm = T)) %>%
    dplyr::select(threshold, confirmation_lens, admin_level, tp_cumu_lb, tp_cumu_median, tp_cumu_ub)
    
df_eff <- tp_ac_allISOs %>%
   filter(year == 2035) %>%
   mutate(efficiency = true_ac_cumu / target_pop_cumu * 1000) %>%
   group_by(threshold, confirmation_lens, admin_level) %>%
   summarize(efficiency_median = median(efficiency, na.rm = TRUE),
             efficiency_lb = quantile(efficiency, 0.025, na.rm = TRUE),
             efficiency_ub = quantile(efficiency, 0.975, na.rm = TRUE)) %>%
   dplyr::select(threshold, confirmation_lens, admin_level, efficiency_lb, efficiency_median, efficiency_ub)


df_ac <- tp_ac_allISOs %>%
    filter(year == 2035) %>%
    group_by(threshold, confirmation_lens, admin_level) %>%
    summarize(ac_median = median(true_ac_cumu, na.rm = TRUE),
              ac_lb = quantile(true_ac_cumu, 0.025, na.rm = TRUE),
              ac_ub = quantile(true_ac_cumu, 0.975, na.rm = TRUE)) %>%
    dplyr::select(threshold, confirmation_lens, admin_level, ac_lb, ac_median, ac_ub)

# join tables 
tp_eff_ac_allISOs_medianCI <- df_tp %>%
    full_join(df_eff, by = c("threshold", "confirmation_lens", "admin_level")) %>%
    full_join(df_ac, by = c("threshold", "confirmation_lens", "admin_level"))
    
write.csv(tp_eff_ac_allISOs_medianCI, paste0(output_final_path, "/intermediate_table/tp_eff_ac_allISOs_medianCI.csv"), row.names = F)


# OUTPUT 5 n_targeted_admins_byISO =======================================================
# get number of targets (by country), aiming to get numbers (across all years) of 1) total targets and 2) unique targets:
for(country in all_countries){
    
    df_tt <- read.csv(paste0(output_final_path, "/target_table/target_table_", country, ".csv" ))
    df_tt <- df_tt %>% mutate(true_averted_cases = no_vaccination_true_case - campaign_default_true_case)
    temp <- df_tt %>% filter(is_target == 1 & true_averted_cases != 0)

    if(which(all_countries == country) == 1){df_target <- temp}
    if(which(all_countries == country) != 1){df_target <- rbind(df_target, temp)}
    
}

df_target <- df_target %>%
        mutate(NAME = case_when(admin_level == "admin1" ~ NAME_1,
                                admin_level == "admin2" ~ NAME_2)) 

# by country
n_total_targets <- df_target %>%
    group_by(ISO, confirmation_lens, admin_level, threshold, run_id) %>%
    summarize(n_total_targets = n()) %>% 
    group_by(ISO, confirmation_lens, admin_level, threshold) %>%
    summarize(n_total_targets_lb = round(quantile(n_total_targets, 0.025, na.rm = T)),
              n_total_targets_median = round(quantile(n_total_targets, 0.5, na.rm = T)),
              n_total_targets_ub = round(quantile(n_total_targets, 0.975, na.rm = T)))

n_unique_targets <- df_target %>%
    group_by(ISO, confirmation_lens, admin_level, threshold, run_id) %>% 
    summarize(n_unique_targets = length(unique(NAME))) %>%
    group_by(ISO, confirmation_lens, admin_level, threshold) %>%
    summarize(n_unique_targets_lb = round(quantile(n_unique_targets, 0.025, na.rm = T)),
              n_unique_targets_median = round(quantile(n_unique_targets, 0.5, na.rm = T)),
              n_unique_targets_ub = round(quantile(n_unique_targets, 0.975, na.rm = T)))

n_targeted_admins_byISO <- 
    full_join(n_total_targets, n_unique_targets,
              by = c("ISO", "confirmation_lens", "admin_level", "threshold"))

write.csv(n_targeted_admins_byISO, "output_final/202302_survms/intermediate_table/n_targeted_admins_byISO.csv", row.names = F)


# OUTPUT 6 n_targeted_admins_allISOs =======================================================
# get number of targets (all countries combined), aiming to get numbers (across all years) of 1) total targets and 2) unique targets:
n_total_targets <- df_target %>%
    group_by(confirmation_lens, admin_level, threshold, run_id) %>%
    summarize(n_total_targets = n()) %>% 
    group_by(confirmation_lens, admin_level, threshold) %>%
    summarize(n_total_targets_lb = round(quantile(n_total_targets, 0.025, na.rm = T)),
              n_total_targets_median = round(quantile(n_total_targets, 0.5, na.rm = T)),
              n_total_targets_ub = round(quantile(n_total_targets, 0.975, na.rm = T)))

n_unique_targets <- df_target %>%
    group_by(confirmation_lens, admin_level, threshold, run_id) %>% 
    summarize(n_unique_targets = length(unique(NAME))) %>%
    group_by(confirmation_lens, admin_level, threshold) %>%
    summarize(n_unique_targets_lb = round(quantile(n_unique_targets, 0.025, na.rm = T)),
              n_unique_targets_median = round(quantile(n_unique_targets, 0.5, na.rm = T)),
              n_unique_targets_ub = round(quantile(n_unique_targets, 0.975, na.rm = T)))

n_targeted_admins_allISOs <- 
    full_join(n_total_targets, n_unique_targets,
              by = c("confirmation_lens", "admin_level", "threshold"))

write.csv(n_targeted_admins_allISOs, "output_final/202302_survms/intermediate_table/n_targeted_admins_allISOs.csv", row.names = F)


# OUTPUT 7: df_small_targets =======================================================
# table of small admins (containing 0 full 5*5 grid) are missed by the fasterize() function, like we discussed before
for(country in all_countries){
    
    df_tt <- read.csv(paste0(output_final_path, "/target_table/target_table_", country, ".csv" ))
    df_tt <- df_tt %>% mutate(true_averted_cases = no_vaccination_true_case - campaign_default_true_case)
    temp <- filter(df_tt, is_target == 1 & true_averted_cases == 0)
    
    if(which(all_countries == country) == 1){df_problem <- temp}
    if(which(all_countries == country) != 1){df_problem <- rbind(df_problem, temp)}

}

write.csv(df_problem, "output_final/202302_survms/intermediate_table/df_small_targets.csv", row.names = F)


# OUTPUT 8: df_small_targets_summary =======================================================
# look into the small admins issue raised by fasterize() function and get a summary table of these admins with info like area (km^2), pop etc. for each admin
if(!file.exists(paste0(output_final_path, "/intermediate_table/df_small_targets.csv"))){
    message("Data file df_small_targets.csv does not exist in the intermediate_table folder.")
}else{
    df_problem <- read.csv(paste0(output_final_path, "/intermediate_table/df_small_targets.csv"))
}

# read in pop and incidence data
pop <- raster(paste0("input_data/worldpop/ppp_2020_5km_Aggregated.tif"))
inc <- raster(paste0("input_data/incidence/afro_2010-2016_lambda_5k_mean.tif"))
case <- pop * inc
ISOs <- unique(df_problem$ISO)
districts <- unique(df_problem$NAME_2)

# get area of the problematic admins
shp2 <- gadm_sf_loadCountries(ISOs, level = 2)$sf
shp2 <- shp2 %>% filter(NAME_2 %in% districts)
shp2$area_km2 <- round(as.numeric(st_area(shp2) / 1000000),1)

# get pop(2020), incidence rate (from mean suspected incidence rate raster of the posterior draws used for baseline) of these admins
pop_admin2 <- exact_extract(pop, shp2, "sum")
shp2$pop <- pop_admin2
case_admin2 <- exact_extract(case, shp2, "sum")
shp2$case <- case_admin2
shp2$incid <- shp2$case / shp2$pop

df_demo <- st_drop_geometry(shp2) %>%
  dplyr::select(ISO, NAME_1, NAME_2, area_km2, pop, case, incid)

# then look at the scenarios these districts were targeted
# (first, note that all these were admin2 targeting)
df_confirmation <- df_problem %>%
  group_by(ISO, NAME_1, NAME_2, confirmation_lens) %>%
  count() %>%
  pivot_wider(names_from = "confirmation_lens", values_from = "n") %>%
  mutate(`global-estimate` = case_when(!is.na(`global-estimate`) ~ "X"),
         `no-estimate` = case_when(!is.na(`no-estimate`) ~ "X"),
         `district-estimate` = case_when(!is.na(`district-estimate`) ~ "X")
         )
  
df_threshold <- df_problem %>%
  group_by(ISO, NAME_1, NAME_2, threshold) %>%
  count() %>%
  pivot_wider(names_from = "threshold", values_from = "n") %>%
  mutate(`1e-04` = case_when(!is.na(`1e-04`) ~ "X"),
         `2e-04` = case_when(!is.na(`2e-04`) ~ "X"),
         `0.001` = case_when(!is.na(`0.001`) ~ "X")
         ) 

df_small_targets_summary <- df_confirmation %>%
  left_join(df_threshold, by = c("ISO", "NAME_1", "NAME_2")) %>%
  left_join(df_demo, by = c("ISO", "NAME_1", "NAME_2")) %>%
  dplyr::select(1:3, 11, 13, 12, 10, 4:9)

write.csv(df_small_targets_summary, "output_final/202302_survms/intermediate_table/df_small_targets_summary.csv", row.names = FALSE)


# OUTPUT 9: n_targeted_runs =======================================================
# number of runs with targeting going on for each scenario by country
for(country in all_countries){
    
    df_tt <- read.csv(paste0(output_final_path, "/target_table/target_table_", country, ".csv" ))
    df_tt <- df_tt %>% mutate(true_averted_cases = no_vaccination_true_case - campaign_default_true_case)
    # remove those small districts that were targeted
    df_tt <- df_tt %>% filter(is_target != 1 | true_averted_cases != 0)
    
    temp_runs <- df_tt %>%
        group_by(ISO, confirmation_lens, threshold, admin_level, run_id) %>%
        summarize(sum_target = sum(is_target)) %>%
        group_by(ISO, confirmation_lens, threshold, admin_level) %>%
        count(sum_target != 0) %>%
        filter(`sum_target != 0` == TRUE) %>%
        rename(n_runs_targeted = n) %>%
        dplyr::select(ISO, threshold, confirmation_lens, admin_level, n_runs_targeted)

    if(which(all_countries == country) == 1){n_targeted_runs <- temp_runs}
    if(which(all_countries == country) != 1){n_targeted_runs <- rbind(n_targeted_runs, temp_runs)}
    
}

write.csv(n_targeted_runs, "output_final/202302_survms/intermediate_table/n_targeted_runs.csv", row.names = F)


# OUTPUT 10: ind_targeted_by_scenario =======================================================
# look at whether a country has any targeting going on by scenario
df_result <- expand_grid(
    threshold = c(1e-04, 2e-04, 1e-03),
    confirmation_lens = c("district-estimate", "global-estimate", "no-estimate"),
    admin_level = c("admin1", "admin2")
)

for(country in all_countries){
    # read in target table 
    df_tt <- read.csv(paste0(output_final_path, "/target_table/target_table_", country, ".csv" ))
    df_tt <- df_tt %>% mutate(true_averted_cases = no_vaccination_true_case - campaign_default_true_case)

    # added on 11/10: remove targets smaller than 5*5 grid (averted case == 0 & is_target == 1)
    df_tt <- df_tt %>% filter(true_averted_cases != 0 | is_target != 1)

    # see if this country has any targeting going on for all the 18 scenarios
    temp <- df_tt %>% filter(is_target == 1) %>% group_by(threshold, confirmation_lens, admin_level) %>% count() %>%
        right_join(df_result, by = c("threshold", "confirmation_lens", "admin_level")) %>%
        mutate(n = ifelse(!is.na(n), 1, 0))
    colnames(temp)[4] <- country

    if(which(all_countries == country) == 1){df_targeted <- temp}
    if(which(all_countries == country) != 1){df_targeted <- df_targeted %>% left_join(temp,  by = c("threshold", "confirmation_lens", "admin_level"))}

}

# if there is no targeting going on in this country for all the simulations (for a given scenario), then noted as 0, otherwise 1.
write.csv(df_targeted, paste0(output_final_path, "/intermediate_table/ind_targeted_by_scenario.csv"), row.names = FALSE)


# OUTPUT 11: ind_targeted_by_run =======================================================
## Then look at for each simulation (run + scenario), whether a country has been targeted
df_result <- expand_grid(
    run_id = 1:200,
    threshold = c(1e-04, 2e-04, 1e-03),
    confirmation_lens = c("district-estimate", "global-estimate", "no-estimate"),
    admin_level = c("admin1", "admin2")
)

for(country in all_countries){
     # read in target table 
    df_tt <- read.csv(paste0(output_final_path, "/target_table/target_table_", country, ".csv" ))
    df_tt <- df_tt %>% mutate(true_averted_cases = no_vaccination_true_case - campaign_default_true_case)

    # added on 11/10: remove targets smaller than 5*5 grid (averted case == 0 & is_target == 1)
    df_tt <- df_tt %>% filter(true_averted_cases != 0 | is_target != 1)

    temp <- df_tt %>% filter(is_target == 1) %>% group_by(run_id, threshold, confirmation_lens, admin_level) %>%
        count() %>%
        full_join(df_result, by = c("run_id", "threshold", "confirmation_lens", "admin_level")) %>%
        mutate(ind_targeted = ifelse(!is.na(n), 1, 0)) %>%
        dplyr::select(-n)
    colnames(temp)[5] <- country
    
    if(which(all_countries == country) == 1){df_targeted_by_run <- temp}
    if(which(all_countries == country) != 1){df_targeted_by_run <- df_targeted_by_run %>% left_join(temp,  by = c("run_id", "threshold", "confirmation_lens", "admin_level"))}
}

write.csv(df_targeted_by_run, paste0(output_final_path, "/intermediate_table/ind_targeted_by_run.csv"), row.names = FALSE)


# OUTPUT 12: n_targeted_admins_ISOs =======================================================
# Number of admins and ISO targeted by scenario 
n_targeted_ISOs <- df_targeted_by_run %>%
  pivot_longer(cols = 5:39, names_to = "ISO", values_to = "ind_targeted") %>%
  group_by(confirmation_lens, threshold, admin_level, run_id) %>%
  summarize(n_targeted_ISO = sum(ind_targeted))  %>%
  group_by(confirmation_lens, threshold, admin_level) %>%
  summarize(n_targeted_ISO_lb = quantile(n_targeted_ISO, 0.025, na.rm = T),
            n_targeted_ISO_median = quantile(n_targeted_ISO, 0.5, na.rm = T),
            n_targeted_ISO_ub = quantile(n_targeted_ISO, 0.975, na.rm = T)) %>%
  mutate(n_targeted_ISO_lb = round(n_targeted_ISO_lb),
         n_targeted_ISO_median = round(n_targeted_ISO_median),
         n_targeted_ISO_ub = round(n_targeted_ISO_ub))

# combine with targeted admin units 
n_targeted_admins_ISOs <- n_targeted_admins_allISOs %>%
    left_join(n_targeted_ISOs, by = c("confirmation_lens", "threshold", "admin_level"))

write.csv(n_targeted_admins_ISOs, "output_final/202302_survms/intermediate_table/n_targeted_admins_ISOs.csv", row.names = FALSE)

