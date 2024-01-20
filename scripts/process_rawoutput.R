# Script to summarize the raw output to get intermediate output datasets
# the intermediate datasets are then to be used in the rmd to make the manuscript figures/tables

# set up =======================================================
r_lib <- Sys.getenv("R_LIBRARY_DIRECTORY", FALSE)
# r_lib <- "/home/hxu70/rlibs/gcm/4.0.2/gcc/9.3.0/"

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
source("packages/ocvImpact/R/utils_process_rawoutput.R")
output_final_path <- "output_final/202110gavi-3"
all_countries <- c("AGO", "BDI", "BEN", "BFA", "CAF", "CIV", "CMR", "COD", "COG", "ETH", "GHA", "GIN", "GNB", "KEN", "LBR", "MDG", "MLI", "MOZ", "MRT", "MWI", "NAM", "NER", "NGA", "RWA", "SEN", "SLE", "SOM", "SSD", "TCD", "TGO", "TZA", "UGA", "ZAF", "ZMB", "ZWE")

# set output paths ===========
int_path <- paste0(output_final_path, "/intermediate_table/")
cost_path <- paste0(output_final_path, "/cost_table/")
tt_path <- paste0(output_final_path, "/target_table/")



# check directory status =======================================================
if(!file.exists(paste0(output_final_path, "/", "target_table"))){
    dir.create(paste0(output_final_path, "/", "target_table"))
}
if(!file.exists(paste0(output_final_path, "/", "intermediate_table"))){
    dir.create(paste0(output_final_path, "/", "intermediate_table"))
}

# OUTPUT 1: target_table =======================================================
# make target_table for each country and store them separately (one country per csv) into "output_final/202110gavi-3/target_table", this table combines all scenarios together for each country
#for(country in all_countries){
#    
#    fn <- paste0("target_table_", country, ".csv")
#    target_table  <- make_target_table(countries = country, thresholds = c("1e-04", "2e-04", "0.001"), output_raw_path = "output_raw/202110gavi-3")
#    write.csv(target_table, paste0(output_final_path, "/target_table/", fn), row.names = F)
#    
#}


# OUTPUT 2: tp_eff_ac_byISO_medianCI table =======================================================
# make table of median and CI of targeted population, OCV efficiency and averted true cases by country
df_result <- expand_grid( # for output 10
    threshold = c(1e-04, 2e-04, 1e-03),
    confirmation_lens = c("district-estimate", "global-estimate", "no-estimate"),
    admin_level = c("admin1", "admin2")
)
df_result_byrun <- expand_grid( # for output 11
    run_id = 1:200,
    threshold = c(1e-04, 2e-04, 1e-03),
    confirmation_lens = c("district-estimate", "global-estimate", "no-estimate"),
    admin_level = c("admin1", "admin2")
)

for(country in all_countries){
    
    # read in target table 
    df_tt <- read.csv(paste0(output_final_path, "/target_table/target_table_", country, ".csv" ))
    df_tt <- df_tt %>% mutate(true_averted_cases = no_vaccination_true_case - campaign_default_true_case)

    # OUTPUT 7:
    temp <- filter(df_tt, is_target == 1 & true_averted_cases == 0)
    if(which(all_countries == country) == 1){df_problem <- temp}
    if(which(all_countries == country) != 1){df_problem <- rbind(df_problem, temp)}

    # added on 11/10/2022: remove targets smaller than 5*5 grid (averted case == 0 & is_target == 1)
    df_tt <- df_tt %>% filter(true_averted_cases != 0 | is_target != 1)
    
    # OUTPUT 2:
    # get targeted pop and efficiency table for this one country
    df_tp_eff <- combine_eff_table_median(countries = country, df_target = df_tt, admin_level = "both")
    if(which(all_countries == country) == 1){df_tp_eff_byISO_medianCI <- df_tp_eff}
    if(which(all_countries == country) != 1){df_tp_eff_byISO_medianCI <- rbind(df_tp_eff_byISO_medianCI, df_tp_eff)}

    # OUTPUT 3:
    if(which(all_countries == country) == 1){df_tt_allISOs <- df_tt}
    if(which(all_countries == country) != 1){df_tt_allISOs <- rbind(df_tt_allISOs, df_tt)}

    # OUTPUT 9
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

    # OUTPUT 10: 
    # see if this country has any targeting going on for all the 18 scenarios
    targeted_scenarios_temp <- df_tt %>% filter(is_target == 1) %>% group_by(threshold, confirmation_lens, admin_level) %>% count() %>%
        right_join(df_result, by = c("threshold", "confirmation_lens", "admin_level")) %>%
        mutate(n = ifelse(!is.na(n), 1, 0))
    colnames(targeted_scenarios_temp)[4] <- country

    if(which(all_countries == country) == 1){df_targeted_scenarios <- targeted_scenarios_temp}
    if(which(all_countries == country) != 1){df_targeted_scenarios <- df_targeted_scenarios %>% left_join(targeted_scenarios_temp,  by = c("threshold", "confirmation_lens", "admin_level"))}

    # OUTPUT 11:
    targeted_runs_temp <- df_tt %>% filter(is_target == 1) %>% group_by(run_id, threshold, confirmation_lens, admin_level) %>%
        count() %>%
        full_join(df_result_byrun, by = c("run_id", "threshold", "confirmation_lens", "admin_level")) %>%
        mutate(ind_targeted = ifelse(!is.na(n), 1, 0)) %>%
        dplyr::select(-n)
    colnames(targeted_runs_temp)[5] <- country
    
    if(which(all_countries == country) == 1){df_targeted_by_run <- targeted_runs_temp}
    if(which(all_countries == country) != 1){df_targeted_by_run <- df_targeted_by_run %>% left_join(targeted_runs_temp,  by = c("run_id", "threshold", "confirmation_lens", "admin_level"))}

}

# OUTPUT 5:
df_target <- df_tt_allISOs

# write output 7
write.csv(df_problem, "output_final/202110gavi-3/intermediate_table/df_small_targets.csv", row.names = F)

# write output 2
if(!file.exists(paste0(output_final_path, "/intermediate_table"))){
    dir.create(paste0(output_final_path, "/intermediate_table"))
    write.csv(df_tp_eff_byISO_medianCI, paste0(output_final_path, "/intermediate_table/tp_eff_ac_byISO_medianCI.csv"), row.names = F)
}else{
    write.csv(df_tp_eff_byISO_medianCI, paste0(output_final_path, "/intermediate_table/tp_eff_ac_byISO_medianCI.csv"), row.names = F)
}

# write output 9
write.csv(n_targeted_runs, "output_final/202110gavi-3/intermediate_table/n_targeted_runs.csv", row.names = F)

# write output 10
write.csv(df_targeted_scenarios, paste0(output_final_path, "/intermediate_table/ind_targeted_by_scenario.csv"), row.names = FALSE)

# write output 11
write.csv(df_targeted_by_run, paste0(output_final_path, "/intermediate_table/ind_targeted_by_run.csv"), row.names = FALSE)


# OUTPUT 3: tp_ac_allISOs =======================================================
# make table of targeted population, efficiency and averted true cases for all countries combined together (summarize across Africa)
# SEE ABOVE (combined processing in one for loop)
#for(country in all_countries){ # stack target tables of all countries
#
#    df_tt <- read.csv(paste0(output_final_path, "/target_table/target_table_", country, ".csv" ))
#    if(which(all_countries == country) == 1){df_tt_allISOs <- df_tt}
#    if(which(all_countries == country) != 1){df_tt_allISOs <- rbind(df_tt_allISOs, df_tt)}
#
#}

# write output 3:
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
# note the target_pop_cumu and true_ac_cumu columns in tp_ac_allISOs.csv are already cumulative by year (rather than value of that year), e.g. if year == 2023, then target_pop_cumu is the target_pop of 2022+2023.



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


# SEE ABOVE (combined processing in one for loop)
# OUTPUT 5 n_targeted_admins_byISO =======================================================
# get number of targets (by country), aiming to get numbers (across all years) of 1) total targets and 2) unique targets:
# for(country in all_countries){
#    
#    df_tt <- read.csv(paste0(output_final_path, "/target_table/target_table_", country, ".csv" ))
#    df_tt <- df_tt %>% mutate(true_averted_cases = no_vaccination_true_case - campaign_default_true_case)
#    temp <- df_tt %>% filter(is_target == 1 & true_averted_cases != 0)
#
#    if(which(all_countries == country) == 1){df_target <- temp}
#    if(which(all_countries == country) != 1){df_target <- rbind(df_target, temp)}
#    
#}

# write ouput 5
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

write.csv(n_targeted_admins_byISO, "output_final/202110gavi-3/intermediate_table/n_targeted_admins_byISO.csv", row.names = F)


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

write.csv(n_targeted_admins_allISOs, "output_final/202110gavi-3/intermediate_table/n_targeted_admins_allISOs.csv", row.names = F)


# OUTPUT 7: df_small_targets =======================================================
# table of small admins (containing 0 full 5*5 grid) are missed by the fasterize() function, like we discussed before
# SEE ABOVE (combined processing in one for loop)
#for(country in all_countries){
#    
#    df_tt <- read.csv(paste0(output_final_path, "/target_table/target_table_", country, ".csv" ))
#    df_tt <- df_tt %>% mutate(true_averted_cases = no_vaccination_true_case - campaign_default_true_case)
#    temp <- filter(df_tt, is_target == 1 & true_averted_cases == 0)
#    
#    if(which(all_countries == country) == 1){df_problem <- temp}
#    if(which(all_countries == country) != 1){df_problem <- rbind(df_problem, temp)}
#
#}
# write.csv(df_problem, "output_final/202110gavi-3/intermediate_table/df_small_targets.csv", row.names = F)


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

write.csv(df_small_targets_summary, "output_final/202110gavi-3/intermediate_table/df_small_targets_summary.csv", row.names = FALSE)


# OUTPUT 9: n_targeted_runs =======================================================
# number of runs with targeting going on for each scenario by country
# SEE ABOVE (combine processing in one for loop)
#for(country in all_countries){
#    
#    df_tt <- read.csv(paste0(output_final_path, "/target_table/target_table_", country, ".csv" ))
#    df_tt <- df_tt %>% mutate(true_averted_cases = no_vaccination_true_case - campaign_default_true_case)
#    # remove those small districts that were targeted
#    df_tt <- df_tt %>% filter(is_target != 1 | true_averted_cases != 0)
#    
#    temp_runs <- df_tt %>%
#        group_by(ISO, confirmation_lens, threshold, admin_level, run_id) %>%
#        summarize(sum_target = sum(is_target)) %>%
#        group_by(ISO, confirmation_lens, threshold, admin_level) %>%
#        count(sum_target != 0) %>%
#        filter(`sum_target != 0` == TRUE) %>%
#        rename(n_runs_targeted = n) %>%
#        dplyr::select(ISO, threshold, confirmation_lens, admin_level, n_runs_targeted)
#
#    if(which(all_countries == country) == 1){n_targeted_runs <- temp_runs}
#    if(which(all_countries == country) != 1){n_targeted_runs <- rbind(n_targeted_runs, temp_runs)}
#    
#}
#write.csv(n_targeted_runs, "output_final/202110gavi-3/intermediate_table/n_targeted_runs.csv", row.names = F)


# OUTPUT 10: ind_targeted_by_scenario =======================================================
# look at whether a country has any targeting going on by scenario
# SEE ABOVE (combine processing code together)
#df_result <- expand_grid(
#    threshold = c(1e-04, 2e-04, 1e-03),
#    confirmation_lens = c("district-estimate", "global-estimate", "no-estimate"),
#    admin_level = c("admin1", "admin2")
#)

#for(country in all_countries){
#    # read in target table 
#    df_tt <- read.csv(paste0(output_final_path, "/target_table/target_table_", country, ".csv" ))
#    df_tt <- df_tt %>% mutate(true_averted_cases = no_vaccination_true_case - campaign_default_true_case)
#
#    # added on 11/10: remove targets smaller than 5*5 grid (averted case == 0 & is_target == 1)
#    df_tt <- df_tt %>% filter(true_averted_cases != 0 | is_target != 1)
#
#    # see if this country has any targeting going on for all the 18 scenarios
#    targeted_scenarios_temp <- df_tt %>% filter(is_target == 1) %>% group_by(threshold, confirmation_lens, admin_level) %>% count() %>%
#        right_join(df_result, by = c("threshold", "confirmation_lens", "admin_level")) %>%
#        mutate(n = ifelse(!is.na(n), 1, 0))
#    colnames(targeted_scenarios_temp)[4] <- country
#
#    if(which(all_countries == country) == 1){df_targeted_scenarios <- targeted_scenarios_temp}
#    if(which(all_countries == country) != 1){df_targeted_scenarios <- df_targeted_scenarios %>% left_join(targeted_scenarios_temp,  by = c("threshold", "confirmation_lens", "admin_level"))}
#
#}

# if there is no targeting going on in this country for all the simulations (for a given scenario), then noted as 0, otherwise 1.
# write.csv(df_targeted_scenarios, paste0(output_final_path, "/intermediate_table/ind_targeted_by_scenario.csv"), row.names = FALSE)


# OUTPUT 11: ind_targeted_by_run =======================================================
# SEE ABOVE(combine processing code in one for loop)
## Then look at for each simulation (run + scenario), whether a country has been targeted
#df_result <- expand_grid(
#    run_id = 1:200,
#    threshold = c(1e-04, 2e-04, 1e-03),
#    confirmation_lens = c("district-estimate", "global-estimate", "no-estimate"),
#    admin_level = c("admin1", "admin2")
#)
#
#for(country in all_countries){
#     # read in target table 
#    df_tt <- read.csv(paste0(output_final_path, "/target_table/target_table_", country, ".csv" ))
#    df_tt <- df_tt %>% mutate(true_averted_cases = no_vaccination_true_case - campaign_default_true_case)
#
#    # added on 11/10: remove targets smaller than 5*5 grid (averted case == 0 & is_target == 1)
#    df_tt <- df_tt %>% filter(true_averted_cases != 0 | is_target != 1)
#
#    targeted_runs_temp <- df_tt %>% filter(is_target == 1) %>% group_by(run_id, threshold, confirmation_lens, admin_level) %>%
#        count() %>%
#        full_join(df_result, by = c("run_id", "threshold", "confirmation_lens", "admin_level")) %>%
#        mutate(ind_targeted = ifelse(!is.na(n), 1, 0)) %>%
#        dplyr::select(-n)
#    colnames(targeted_runs_temp)[5] <- country
#    
#    if(which(all_countries == country) == 1){df_targeted_by_run <- targeted_runs_temp}
#    if(which(all_countries == country) != 1){df_targeted_by_run <- df_targeted_by_run %>% left_join(targeted_runs_temp,  by = c("run_id", "threshold", "confirmation_lens", "admin_level"))}
#}
#
#write.csv(df_targeted_by_run, paste0(output_final_path, "/intermediate_table/ind_targeted_by_run.csv"), row.names = FALSE)


# OUTPUT 12: n_targeted_admins_ISOs =======================================================
# Number of admins and ISO targeted by scenario (as well as the number of countries targeted in each scenario)
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

write.csv(n_targeted_admins_ISOs, "output_final/202110gavi-3/intermediate_table/n_targeted_admins_ISOs.csv", row.names = FALSE)



# OUTPUT 24: mean_ir_targeted
# mean true baseline IR among admins that were vaccinated at least once in that scenario
for(country in all_countries){
    
    # read in target table 
    message(paste0("For making dataframe of mean IR, reading the target table of ", country))
    df_tt <- read.csv(paste0(output_final_path, "/target_table/target_table_", country, ".csv" ))
    df_tt <- df_tt %>% mutate(true_averted_cases = no_vaccination_true_case - campaign_default_true_case)

    # added on 11/10/2022: remove targets smaller than 5*5 grid (averted case == 0 & is_target == 1)
    df_tt <- df_tt %>% filter(true_averted_cases != 0 | is_target != 1)

    # admins targeted in each scenario
    unique_targeted_admins_temp <- df_tt %>% 
        mutate(NAME = ifelse(admin_level == "admin1", NAME_1, NAME_2)) %>%
        filter(is_target == 1) %>%
        group_by(ISO, confirmation_lens, admin_level, threshold) %>% # targeted admins in all runs
        summarize(targeted_admins = paste0(unique(NAME), collapse = ", ")) %>%
        mutate(mean_true_ir = NA, mean_confirmed_ir = NA)

    if(nrow(unique_targeted_admins_temp) != 0){

        df_base_inc <- df_tt %>% 
        filter(year == 2022 & run_id == 1) %>%
        mutate(NAME = ifelse(admin_level == "admin1", NAME_1, NAME_2)) %>%
        dplyr::select(ISO, admin_level, NAME, confirmation_lens, threshold, true_incidence_rate, confirmed_incidence_rate)


        for(i in 1:nrow(unique_targeted_admins_temp)){
        
            admin_level_temp <- unique_targeted_admins_temp$admin_level[i]
            confirmation_lens_temp <- unique_targeted_admins_temp$confirmation_lens[i]
            threshold_temp <- unique_targeted_admins_temp$threshold[i]
        
            targeted_admins <- unlist(strsplit(unique_targeted_admins_temp$targeted_admins[i], ", "))
            mean_ir_temp <- df_base_inc %>%
                filter(confirmation_lens == confirmation_lens_temp, admin_level == admin_level_temp, threshold == threshold_temp) %>%
                filter(NAME %in% targeted_admins) %>%
                group_by(admin_level, confirmation_lens, threshold) %>%
                summarize(mean_true_ir = mean(true_incidence_rate, na.rm = T),
                        mean_confirmed_ir = mean(confirmed_incidence_rate, na.rm = T))
        
            # fill in the table 
            unique_targeted_admins_temp$mean_true_ir[i] <- mean_ir_temp$mean_true_ir 
            unique_targeted_admins_temp$mean_confirmed_ir[i] <- mean_ir_temp$mean_confirmed_ir
        }

        # unique_targeted_admins
        if(which(all_countries == country) == 1){unique_targeted_admins <- unique_targeted_admins_temp}
        if(which(all_countries == country) != 1){unique_targeted_admins <- rbind(unique_targeted_admins, unique_targeted_admins_temp)}
    
    }
    
}
write.csv(unique_targeted_admins, "output_final/202110gavi-3/intermediate_table/mean_ir_targeted.csv", row.names = FALSE)



### OUTPUT 25: ir_eff_byISO.csv
# Average district-level baseline observed IR ~ efficiency (for all districts targeted at least once)
if(!file.exists(paste0(output_final_path, "/intermediate_table/mean_ir_targeted.csv"))){
    message("Data file mean_ir_targeted.csv does not exist in the intermediate_table folder.")
}else{
    df_targeted_admins <- read.csv(paste0(output_final_path, "/intermediate_table/mean_ir_targeted.csv"))
}

# make a folder for ir_eff tables
if(!file.exists(paste0(output_final_path, "/", "ir_eff_table"))){
    dir.create(paste0(output_final_path, "/", "ir_eff_table"))
}


# remove SEN and ZAF (without targeting in all scenarios)
all_countries <- c("AGO", "BDI", "BEN", "BFA", "CAF", "CIV", "CMR", "COD", "COG", "ETH", "GHA", "GIN", "GNB", "KEN", "LBR", "MDG", "MLI", "MOZ", "MRT", "MWI", "NAM", "NER", "NGA", "RWA", "SLE", "SOM", "SSD", "TCD", "TGO", "TZA", "UGA", "ZMB", "ZWE")


for(country in all_countries){
    
    df_targeted_admins_temp <- df_targeted_admins %>% filter(country == ISO)
    
    # read in target table 
    message(paste0("Reading the target table of ", country))
    df_tt <- read.csv(paste0(output_final_path, "/target_table/target_table_", country, ".csv" ))
    df_tt <- df_tt %>% mutate(true_averted_cases = no_vaccination_true_case - campaign_default_true_case)

    # added on 11/10/2022: remove targets smaller than 5*5 grid (averted case == 0 & is_target == 1)
    df_tt <- df_tt %>% filter(true_averted_cases != 0 | is_target != 1)

    df_tt <- df_tt %>%
        mutate(NAME = ifelse(admin_level == "admin1", NAME_1, NAME_2))
    
    for(i in 1:nrow(df_targeted_admins_temp)){
        
        message(paste0("Reading in line ", i, " of ", country, " targeted admin table."))
        admins <- c(unlist(strsplit(df_targeted_admins_temp$targeted_admins[i], ", ")))
        n_admins <- length(admins)
        
        threshold_temp <- df_targeted_admins_temp$threshold[i]
        admin_level_temp <- df_targeted_admins_temp$admin_level[i]
        confirmation_lens_temp <- df_targeted_admins_temp$confirmation_lens[i]
        
        df_tt_temp <- df_tt %>%
            filter(threshold == threshold_temp & admin_level == admin_level_temp & confirmation_lens == confirmation_lens_temp) %>%
            filter(NAME %in% admins)

        # baseline observed and clinical incidence rate (average across runs) for each admin
        df_inc <- df_tt_temp %>%
            filter(year == 2022) %>%
            group_by(NAME) %>%
            summarize(mean_observed_ir = mean(confirmed_incidence_rate, na.rm = T),
                      mean_clinical_ir = mean(confirmed_incidence_rate / confirmation_rate, na.rm = T))
        
        # averted cases and FVP for each admin
        df_fvp_ac_eff <- df_tt_temp %>%
            group_by(run_id, NAME) %>% 
            mutate(fvp_cumu = cumsum(actual_fvp), ac_cumu = cumsum(true_averted_cases)) %>%
            filter(year == 2035) %>%
            mutate(eff = ac_cumu / fvp_cumu * 1000) %>%
            filter(fvp_cumu != 0) %>% # filter out eff == Inf 
            group_by(NAME) %>% 
            summarize(eff = median(eff, na.rm = TRUE)) %>%
            as.data.frame()
        
        # combine
        df_inc_eff_temp <- df_inc %>% left_join(df_fvp_ac_eff, by = "NAME") %>%
            mutate(ISO = country, threshold = threshold_temp, confirmation_lens = confirmation_lens_temp, admin_level = admin_level_temp) %>%
            dplyr::select(ISO, threshold, confirmation_lens, admin_level, NAME, mean_observed_ir, mean_clinical_ir, eff)
        
        # rbind for one country
        if(i == 1){df_inc_eff_oneISO <- df_inc_eff_temp}
        if(i != 1){df_inc_eff_oneISO <- rbind(df_inc_eff_oneISO, df_inc_eff_temp)}
    }
    # write data for one country
    message(paste0("Writing ir_eff table for ", "country."))
    write.csv(df_inc_eff_oneISO, paste0(output_final_path, "/ir_eff_table/ir_eff_", country, ".csv"), row.names = FALSE)
    
    message(paste0("Rbinding the ir_eff table of ", country, " to byISO table."))
    if(which(all_countries == country) == 1){ir_eff_byISO <- df_inc_eff_oneISO}
    if(which(all_countries == country) != 1){ir_eff_byISO <- rbind(ir_eff_byISO, df_inc_eff_oneISO)}
    
    # rm
    rm(df_tt, df_tt_temp, df_inc, df_fvp_ac_eff, df_inc_eff_temp, df_inc_eff_oneISO)

}
write.csv(ir_eff_byISO, paste0(int_path, "ir_eff_byISO.csv"), row.names = F)



## OUTPUT 26: py_vaxed_byISO
# Number and proportion of population-year getting vaccinated among districts where true IR exceeds the targeting threshold
if(!file.exists(paste0(output_final_path, "/intermediate_table/", "py_vaxed_table"))){
    dir.create(paste0(output_final_path, "/intermediate_table/", "py_vaxed_table"))
}

for(country in all_countries){
    
    # read in target table 
    message(paste0("For making py_vaxed_byISO, reading the target table of ", country))
    df_tt <- read.csv(paste0(output_final_path, "/target_table/target_table_", country, ".csv" ))
    df_tt <- df_tt %>% mutate(true_averted_cases = no_vaccination_true_case - campaign_default_true_case)

    # added on 11/10/2022: remove targets smaller than 5*5 grid (averted case == 0 & is_target == 1)
    df_tt <- df_tt %>% filter(true_averted_cases != 0 | is_target != 1)
    
    df_py_vaccinated <- df_tt %>% 
        filter(is_target == 1 & true_incidence_rate > threshold) %>%
        group_by(confirmation_lens, admin_level, threshold, run_id) %>%
        summarize(py_vaccinated = sum(actual_fvp))
    
    df_prop_py_vaccinated <- df_tt %>%
        filter(true_incidence_rate > threshold) %>%
        group_by(confirmation_lens, admin_level, threshold, run_id) %>%
        summarize(py = sum(pop_model))
    rm(df_tt)

    # join 
    py_vaccinated_temp <- df_py_vaccinated %>%
        left_join(df_prop_py_vaccinated, by = c("confirmation_lens", "admin_level", "threshold", "run_id")) %>%
        mutate(prop_py_vaccinated = py_vaccinated / py) %>%
        group_by(confirmation_lens, admin_level, threshold) %>%
        summarize(py_vaxed_lb = quantile(py_vaccinated, 0.025, na.rm = T),
                  py_vaxed_median = quantile(py_vaccinated, 0.5, na.rm = T),
                  py_vaxed_ub = quantile(py_vaccinated, 0.975, na.rm = T),
                  prop_py_vaxed_lb = quantile(prop_py_vaccinated, 0.025, na.rm = T),
                  prop_py_vaxed_median = quantile(prop_py_vaccinated, 0.5, na.rm = T),
                  prop_py_vaxed_ub = quantile(prop_py_vaccinated, 0.975, na.rm = T)) %>%
        mutate(ISO = country)
    rm(df_py_vaccinated, df_prop_py_vaccinated)
    
    # save the df for a single country
    write.csv(py_vaccinated_temp, paste0(output_final_path, "/intermediate_table/py_vaxed_table/py_vaxed_", country, ".csv"), row.names = F)

    message("Rbinding py_vaccinated_temp of ", country, " to existing py_vaxed_byISO.")
    if(which(all_countries == country) == 1){py_vaxed_byISO <- py_vaccinated_temp}
    if(which(all_countries == country) != 1){py_vaxed_byISO <- rbind(py_vaxed_byISO, py_vaccinated_temp)}
}
write.csv(py_vaxed_byISO, paste0(output_final_path, "/intermediate_table/py_vaxed_byISO.csv"), row.names = F)


## OUTPUT 27: py_vaxed_allISOs
# Number and proportion of population-year getting vaccinated among districts where true IR exceeds the targeting threshold (summarized across all countries)
for(country in all_countries){
    
    # read in target table 
    message(paste0("For making py_vaxed_allISOs, reading the target table of ", country))
    df_tt <- read.csv(paste0(output_final_path, "/target_table/target_table_", country, ".csv" ))
    df_tt <- df_tt %>% mutate(true_averted_cases = no_vaccination_true_case - campaign_default_true_case)

    # added on 11/10/2022: remove targets smaller than 5*5 grid (averted case == 0 & is_target == 1)
    df_tt <- df_tt %>% filter(true_averted_cases != 0 | is_target != 1)

    # summarize the py_vaxed for each run & scenario
    df_py_vaxed <- df_tt %>%
        filter(is_target == 1 & true_incidence_rate > threshold) %>%
        group_by(confirmation_lens, threshold, admin_level, run_id) %>%
        summarize(py_vaxed = sum(actual_fvp))
    
    df_py <- df_tt %>%
        filter(true_incidence_rate > threshold) %>%
        group_by(confirmation_lens, threshold, admin_level, run_id) %>%
        summarize(py = sum(pop_model))
    rm(df_tt)
    
    py_vaxed_temp <- df_py_vaxed %>%
        left_join(df_py, by = c("confirmation_lens", "threshold", "admin_level", "run_id")) %>%
        mutate(ISO = country)
    write.csv(py_vaxed_temp, paste0(output_final_path, "/intermediate_table/py_vaxed_table/py_vaxed_allISOs_", country, ".csv"), row.names = F)
    
    message("Rbinding py_vaxed_temp of ", country, " to existing py_vaxed_allISOs.")
    if(which(all_countries == country) == 1){py_vaxed_allISOs_temp <- py_vaxed_temp}
    if(which(all_countries == country) != 1){py_vaxed_allISOs_temp <- rbind(py_vaxed_allISOs_temp, py_vaxed_temp)}
    rm(py_vaxed_temp)
}
py_vaxed_allISOs <- py_vaxed_allISOs_temp %>%
    group_by(confirmation_lens, threshold, admin_level, run_id) %>%
    summarize(py_vaxed = sum(py_vaxed),
              py = sum(py)) %>%
    mutate(prop_py_vaxed = py_vaxed / py) %>%
    group_by(confirmation_lens, threshold, admin_level) %>%
    summarize(py_vaxed_lb = quantile(py_vaxed, 0.025, na.rm = T),
                  py_vaxed_median = quantile(py_vaxed, 0.5, na.rm = T),
                  py_vaxed_ub = quantile(py_vaxed, 0.975, na.rm = T),
                  prop_py_vaxed_lb = quantile(prop_py_vaxed, 0.025, na.rm = T),
                  prop_py_vaxed_median = quantile(prop_py_vaxed, 0.5, na.rm = T),
                  prop_py_vaxed_ub = quantile(prop_py_vaxed, 0.975, na.rm = T))
write.csv(py_vaxed_allISOs, paste0(output_final_path, "/intermediate_table/py_vaxed_allISOs.csv"), row.names = F)



## OUTPUT 32: n_prop_ac_allISOs_medianCI
# number of averted true cases, and proportion of averted cases among true cases (summarize across countries)
for(country in all_countries){
        
    # read in target table 
    message(paste0("For making n/prop or averted cases, reading the target table of ", country))
    df_tt <- read.csv(paste0(output_final_path, "/target_table/target_table_", country, ".csv" ))
    df_tt <- df_tt %>% mutate(true_averted_cases = no_vaccination_true_case - campaign_default_true_case)

    # added on 11/10/2022: remove targets smaller than 5*5 grid (averted case == 0 & is_target == 1)
    df_tt <- df_tt %>% filter(true_averted_cases != 0 | is_target != 1)

    # calculate number of true cases in no-vaccination scenario
    df_true_cases_temp <- df_tt %>%
        filter(confirmation_lens == "district-estimate") %>%
        group_by(confirmation_lens, threshold, admin_level, run_id) %>%
        summarize(cumu_true_cases = sum(no_vaccination_true_case))
    write.csv(df_true_cases_temp, paste0(int_path, "true_cases_table/true_cases_", country, ".csv"), row.names = F)
    rm(df_tt)
    
    message(paste0("Rbinding the df_true_cases_temp table of ", country, " to df_true_cases table."))
    if(which(all_countries == country) == 1){df_true_cases <- df_true_cases_temp}
    if(which(all_countries == country) != 1){df_true_cases <- rbind(df_true_cases, df_true_cases_temp)}
    rm(df_true_cases_temp)
}
true_cases_allISOs <- df_true_cases %>%
    group_by(confirmation_lens, threshold, admin_level, run_id) %>%
    summarize(true_cases = sum(cumu_true_cases))

write.csv(true_cases_allISOs, paste0(int_path, "true_cases_allISOs.csv"), row.names = F)

true_cases_allISOs <- read.csv(paste0(int_path, "true_cases_allISOs.csv"))
tp_ac_allISOs <- read.csv(paste0(int_path, "tp_ac_allISOs.csv"))

n_prop_ac_allISOs_medianCI <- tp_ac_allISOs %>%
    filter(year == 2035) %>%
    left_join(
        true_cases_allISOs %>% filter(threshold == 1e-04) %>% dplyr::select(-confirmation_lens, -threshold)
    ) %>%
    mutate(prop_ac = true_ac_cumu / true_cases) %>%
    group_by(confirmation_lens, threshold, admin_level) %>%
    summarize(true_ac_lb = quantile(true_ac_cumu, 0.025, na.rm = T),
              true_ac_median = quantile(true_ac_cumu, 0.5, na.rm = T),
              true_ac_ub = quantile(true_ac_cumu, 0.975, na.rm = T),
              prop_ac_lb = quantile(prop_ac, 0.025, na.rm = T),
              prop_ac_median = quantile(prop_ac, 0.5, na.rm = T),
              prop_ac_ub = quantile(prop_ac, 0.975, na.rm = T))
write.csv(n_prop_ac_allISOs_medianCI, paste0(int_path, "n_prop_ac_allISOs_medianCI.csv"), row.names = F)


## OUTPUT 33 prop_ac_byISO_medianCI
# number of averted true cases, and proportion of averted cases among true cases (by country)
true_cases_path <- paste0(int_path, 'true_cases_table')
prop_ac_path <- paste0(int_path, "prop_ac_table")
for(country in all_countries){
        
    # read in target table 
    message(paste0("For making n/prop or averted cases, reading the target table of ", country))
    df_tt <- read.csv(paste0(output_final_path, "/target_table/target_table_", country, ".csv" ))
    df_tt <- df_tt %>% mutate(true_averted_cases = no_vaccination_true_case - campaign_default_true_case)

    # added on 11/10/2022: remove targets smaller than 5*5 grid (averted case == 0 & is_target == 1)
    df_tt <- df_tt %>% filter(true_averted_cases != 0 | is_target != 1)

    # read in true cases in no-vax scenario
    df_true_cases <- read.csv(paste0(true_cases_path, "/true_cases_", country, ".csv")) %>% mutate(ISO = country)
    
    # get averted cases of the country
    df_ac <- df_tt %>%
        group_by(confirmation_lens, admin_level, threshold, run_id) %>%
        summarize(ac = sum(true_averted_cases))
    rm(df_tt)

    # merge 
    df_prop_ac_temp <- df_ac %>% 
        left_join(df_true_cases %>% dplyr::select(-confirmation_lens), by = c("threshold", "admin_level", "run_id")) %>%
        mutate(prop_ac = ac / cumu_true_cases) %>%
        filter(prop_ac > 0) %>%
        group_by(ISO, confirmation_lens, threshold, admin_level) %>%
        summarize(prop_ac_lb = quantile(prop_ac, 0.025, na.rm = T),
                  prop_ac_median = quantile(prop_ac, 0.5, na.rm = T),
                  prop_ac_ub = quantile(prop_ac, 0.975, na.rm = T))
    rm(df_ac, df_true_cases)
    write.csv(df_prop_ac_temp, paste0(prop_ac_path, "/prop_ac_", country, ".csv"), row.names = F)

    message(paste0("Rbinding the df_prop_ac_temp table of ", country, " to prop_ac_byISO_medianCI table."))
    if(which(all_countries == country) == 1){prop_ac_byISO_medianCI <- df_prop_ac_temp}
    if(which(all_countries == country) != 1){prop_ac_byISO_medianCI <- rbind(prop_ac_byISO_medianCI, df_prop_ac_temp)}

}
write.csv(prop_ac_byISO_medianCI, paste0(int_path, "prop_ac_byISO_medianCI.csv"), row.names = F)




## OUTPUT 38: prop_vaxed_in_highIR_allISOs_medianCI
# Proportion of population-year getting vaccinated among districts where true IR exceeds the targeting threshold (summarized across all countries)
for(country in all_countries){
    
    # read in target table 
    message(paste0("For making py_vaxed_allISOs, reading the target table of ", country))
    df_tt <- read.csv(paste0(output_final_path, "/target_table/target_table_", country, ".csv" ))
    df_tt <- df_tt %>% mutate(true_averted_cases = no_vaccination_true_case - campaign_default_true_case)

    # added on 11/10/2022: remove targets smaller than 5*5 grid (averted case == 0 & is_target == 1)
    df_tt <- df_tt %>% filter(true_averted_cases != 0 | is_target != 1)

    # summarize targeted pop living in high burden areas for each run & scenario
    df_py_vaxed_in_highIR <- df_tt %>%
        filter(is_target == 1 & true_incidence_rate > threshold) %>%
        group_by(ISO, confirmation_lens, threshold, admin_level, run_id) %>%
        summarize(py_vaxed_in_highIR = sum(actual_fvp))
    rm(df_tt)
    write.csv(df_py_vaxed_in_highIR, paste0(output_final_path, "/intermediate_table/py_vaxed_in_highIR_table/py_vaxed_in_highIR_", country, ".csv"), row.names = F)
    
    message("Rbinding py_vaxed_in_highIR of ", country, " to existing py_vaxed_in_highIR_byISO.")
    if(which(all_countries == country) == 1){df_py_vaxed_in_highIR_byISO <- df_py_vaxed_in_highIR}
    if(which(all_countries == country) != 1){df_py_vaxed_in_highIR_byISO <- rbind(df_py_vaxed_in_highIR_byISO, df_py_vaxed_in_highIR)}
    rm(df_py_vaxed_in_highIR)
}
py_vaxed_in_highIR_allISOs <- df_py_vaxed_in_highIR_byISO %>%
    group_by(confirmation_lens, threshold, admin_level, run_id) %>%
    summarize(py_vaxed_in_highIR = sum(py_vaxed_in_highIR)) 
    
py_vaxed_allISOs <- read.csv(paste0(int_path, "tp_ac_allISOs.csv"))

n_vaxed_in_highIR_allISOs <- py_vaxed_in_highIR_allISOs %>%
    left_join(py_vaxed_allISOs %>% filter(year == 2035) %>% dplyr::select(-year), 
              by = c("run_id", "threshold", "confirmation_lens", "admin_level"))
write.csv(n_vaxed_in_highIR_allISOs, paste0(int_path, "n_vaxed_in_highIR_allISOs.csv"), row.names = F)

prop_vaxed_in_highIR_allISOs_medianCI <- n_vaxed_in_highIR_allISOs %>%
    mutate(prop_vaxed_in_highIR = py_vaxed_in_highIR / target_pop_cumu) %>%
    group_by(confirmation_lens, threshold, admin_level) %>%
    summarize(prop_vaxed_in_highIR_lb = quantile(prop_vaxed_in_highIR, 0.025, na.rm = T),
              prop_vaxed_in_highIR_median = quantile(prop_vaxed_in_highIR, 0.5, na.rm = T),
              prop_vaxed_in_highIR_ub = quantile(prop_vaxed_in_highIR, 0.975, na.rm = T))

write.csv(prop_vaxed_in_highIR_allISOs_medianCI, paste0(output_final_path, "/intermediate_table/prop_vaxed_in_highIR_allISOs_medianCI.csv"), row.names = F)


## OUTPUT 39:  prop_vaxed_in_highIR_byISO_medianCI (updated 1.17.2024)
# Proportion of population-year getting vaccinated among districts where true IR exceeds the targeting threshold
for(country in all_countries){
    
    if(!country %in% c("SEN", "ZAF")){
        
        # read in target table 
        message(paste0("For making prop_vaxed_in_highIR_allISOs_medianCI (output 39), reading the target table of ", country))
        df_tt <- read.csv(paste0(output_final_path, "/target_table/target_table_", country, ".csv" ))
        df_tt <- df_tt %>% mutate(true_averted_cases = no_vaccination_true_case - campaign_default_true_case)

        # added on 11/10/2022: remove targets smaller than 5*5 grid (averted case == 0 & is_target == 1)
        df_tt <- df_tt %>% filter(true_averted_cases != 0 | is_target != 1)
        
        # calculate targeted pop living in high burden areas for each run & scenario
        df_py_vaxed_in_highIR <- df_tt %>%
            filter(is_target == 1 & true_incidence_rate > threshold) %>%
            group_by(ISO, confirmation_lens, threshold, admin_level, run_id) %>%
            mutate(py_vaxed_in_highIR = cumsum(actual_fvp)) %>%
            filter(year == 2030) %>%
            dplyr::select(ISO, confirmation_lens, threshold, admin_level, run_id, py_vaxed_in_highIR) %>%
            group_by(ISO, confirmation_lens, threshold, admin_level, run_id) %>%
            summarize(py_vaxed_in_highIR = sum(py_vaxed_in_highIR))

        # calculate targeted population
        df_py_vaxed_all <- df_tt %>%
            filter(is_target == 1) %>%
            group_by(ISO, confirmation_lens, threshold, admin_level, run_id) %>%
            mutate(py_vaxed = cumsum(actual_fvp)) %>%
            filter(year == 2030) %>%
            dplyr::select(ISO, confirmation_lens, threshold, admin_level, run_id, py_vaxed) %>%
            group_by(ISO, confirmation_lens, threshold, admin_level, run_id) %>%
            summarize(py_vaxed = sum(py_vaxed))

        # join together 
        prop_vaxed_in_highIR_oneISO <- df_py_vaxed_all %>%
            full_join(df_py_vaxed_in_highIR, by = c("ISO", "confirmation_lens", "threshold", "admin_level", "run_id")) %>%
            mutate(prop_vaxed_in_highIR = py_vaxed_in_highIR / py_vaxed) %>%
            group_by(ISO, confirmation_lens, threshold, admin_level) %>%
            summarize(prop_vaxed_in_highIR_median = median(prop_vaxed_in_highIR, na.rm = T),
                      prop_vaxed_in_highIR_lb = quantile(prop_vaxed_in_highIR, 0.025, na.rm = T),
                      prop_vaxed_in_highIR_ub = quantile(prop_vaxed_in_highIR, 0.975, na.rm = T)) %>%
            mutate(prop_vaxed_in_highIR_median = ifelse(is.na(prop_vaxed_in_highIR_median), 0, prop_vaxed_in_highIR_median),
                   prop_vaxed_in_highIR_lb = ifelse(is.na(prop_vaxed_in_highIR_lb), 0, prop_vaxed_in_highIR_lb),
                   prop_vaxed_in_highIR_ub = ifelse(is.na(prop_vaxed_in_highIR_ub), 0, prop_vaxed_in_highIR_ub))
        # for scenarios not presented in this above table: no targeting for this scenario
        # for scenarios with NA: there is targeting, but there is not a single run that has true incidence above threshold, thus can be replaced by 0
        write.csv(prop_vaxed_in_highIR_oneISO, paste0(output_final_path, "/intermediate_table/prop_vaxed_in_highIR_table/prop_vaxed_in_highIR_", country, ".csv"), row.names = F)   

        if(which(all_countries == country) == 1){prop_vaxed_in_highIR_byISO <- prop_vaxed_in_highIR_oneISO}
        if(which(all_countries == country) != 1){prop_vaxed_in_highIR_byISO <- rbind(prop_vaxed_in_highIR_byISO, prop_vaxed_in_highIR_oneISO)}
    }

}
message("Writing output 39")
write.csv(prop_vaxed_in_highIR_byISO, paste0(int_path, "prop_vaxed_in_highIR_byISO.csv"), row.names = F)


### OUTPUT 40: targeted_baseline_true_ir_byISO.csv (unused)
# Average district-level baseline true IR (for all districts targeted at least once)
if(!file.exists(paste0(output_final_path, "/intermediate_table/mean_ir_targeted.csv"))){
    message("Data file mean_ir_targeted.csv does not exist in the intermediate_table folder.")
}else{
    df_targeted_admins <- read.csv(paste0(output_final_path, "/intermediate_table/mean_ir_targeted.csv"))
}

# make a folder for true_ir tables
if(!file.exists(paste0(output_final_path, "/", "true_ir_table"))){
    dir.create(paste0(output_final_path, "/", "true_ir_table"))
}


# remove SEN and ZAF (without targeting in all scenarios)
all_countries <- c("AGO", "BDI", "BEN", "BFA", "CAF", "CIV", "CMR", "COD", "COG", "ETH", "GHA", "GIN", "GNB", "KEN", "LBR", "MDG", "MLI", "MOZ", "MRT", "MWI", "NAM", "NER", "NGA", "RWA", "SLE", "SOM", "SSD", "TCD", "TGO", "TZA", "UGA", "ZMB", "ZWE")


for(country in all_countries){
    
    df_targeted_admins_temp <- df_targeted_admins %>% filter(country == ISO)
    
    # read in target table 
    message(paste0("Reading the target table of ", country))
    df_tt <- read.csv(paste0(output_final_path, "/target_table/target_table_", country, ".csv" ))
    df_tt <- df_tt %>% mutate(true_averted_cases = no_vaccination_true_case - campaign_default_true_case)

    # added on 11/10/2022: remove targets smaller than 5*5 grid (averted case == 0 & is_target == 1)
    df_tt <- df_tt %>% filter(true_averted_cases != 0 | is_target != 1)

    df_tt <- df_tt %>%
        mutate(NAME = ifelse(admin_level == "admin1", NAME_1, NAME_2))
    
    for(i in 1:nrow(df_targeted_admins_temp)){
        
        message(paste0("Reading in line ", i, " of ", country, " targeted admin table."))
        admins <- c(unlist(strsplit(df_targeted_admins_temp$targeted_admins[i], ", ")))
        n_admins <- length(admins)
        
        threshold_temp <- df_targeted_admins_temp$threshold[i]
        admin_level_temp <- df_targeted_admins_temp$admin_level[i]
        confirmation_lens_temp <- df_targeted_admins_temp$confirmation_lens[i]
        
        df_tt_temp <- df_tt %>%
            filter(threshold == threshold_temp & admin_level == admin_level_temp & confirmation_lens == confirmation_lens_temp) %>%
            filter(NAME %in% admins)

        # baseline observed and clinical incidence rate (average across runs) for each admin
        df_inc <- df_tt_temp %>%
            filter(year == 2022) %>%
            group_by(ISO, threshold, confirmation_lens, admin_level, NAME) %>%
            summarize(true_baseline_ir = mean(true_incidence_rate, na.rm = T)) 
        
        # rbind for one country
        if(i == 1){df_inc_oneISO <- df_inc}
        if(i != 1){df_inc_oneISO <- rbind(df_inc_oneISO, df_inc)}
    }
    # write data for one country
    message(paste0("Writing true IR table for ", country))
    write.csv(df_inc_oneISO, paste0(output_final_path, "/true_ir_table/true_ir_", country, ".csv"), row.names = FALSE)
    
    message(paste0("Rbinding the true IR table of ", country, " to byISO table."))
    if(which(all_countries == country) == 1){targeted_baseline_true_ir_byISO <- df_inc_oneISO}
    if(which(all_countries == country) != 1){targeted_baseline_true_ir_byISO <- rbind(targeted_baseline_true_ir_byISO, df_inc_oneISO)}
    
    # rm
    rm(df_tt, df_tt_temp, df_inc, df_inc_oneISO)

}
write.csv(targeted_baseline_true_ir_byISO, paste0(int_path, "targeted_baseline_true_ir_byISO.csv"), row.names = F)



### OUTPUT 40: targeted_baseline_true_ir_byISO.csv
# Average district-level baseline true IR (for all districts targeted at least once)
if(!file.exists(paste0(output_final_path, "/intermediate_table/mean_ir_targeted.csv"))){
    message("Data file mean_ir_targeted.csv does not exist in the intermediate_table folder.")
}else{
    df_targeted_admins <- read.csv(paste0(output_final_path, "/intermediate_table/mean_ir_targeted.csv"))
}

# make a folder for true_ir tables
if(!file.exists(paste0(output_final_path, "/", "true_ir_table"))){
    dir.create(paste0(output_final_path, "/", "true_ir_table"))
}


# remove SEN and ZAF (without targeting in all scenarios)
all_countries <- c("AGO", "BDI", "BEN", "BFA", "CAF", "CIV", "CMR", "COD", "COG", "ETH", "GHA", "GIN", "GNB", "KEN", "LBR", "MDG", "MLI", "MOZ", "MRT", "MWI", "NAM", "NER", "NGA", "RWA", "SLE", "SOM", "SSD", "TCD", "TGO", "TZA", "UGA", "ZMB", "ZWE")


for(country in all_countries){
        
    # read in target table 
    message(paste0("Reading the target table of ", country))
    df_tt <- read.csv(paste0(output_final_path, "/target_table/target_table_", country, ".csv" ))
    df_tt <- df_tt %>% mutate(true_averted_cases = no_vaccination_true_case - campaign_default_true_case)

    # added on 11/10/2022: remove targets smaller than 5*5 grid (averted case == 0 & is_target == 1)
    df_tt <- df_tt %>% filter(true_averted_cases != 0 | is_target != 1)
    
    # filter out data rows targeted
    df_inc <- df_tt %>% 
        filter(is_target == 1 & admin_level == "admin2") %>%
        filter(confirmation_lens == "district-estimate" | confirmation_lens == "no-estimate") %>% 
        dplyr::select(ISO, admin_level,NAME_2, confirmation_lens, threshold, year, run_id, true_incidence_rate)

    
    # write data for one country
    message(paste0("Writing true IR table for ", country))
    write.csv(df_inc, paste0(output_final_path, "/true_ir_table/true_ir_", country, ".csv"), row.names = FALSE)
    
    # rm
    rm(df_tt, df_inc)

}
write.csv(targeted_baseline_true_ir_byISO, paste0(int_path, "targeted_baseline_true_ir_byISO.csv"), row.names = F)
