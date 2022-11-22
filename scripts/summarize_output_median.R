# script to summarize the raw output

# set up
library(dplyr)
library(stringr)
library(tidyr)
source("packages/ocvImpact/R/utils_sum_output_median.R")
output_final_path <- "output_final/202110gavi-3"

# 10/12/22
# make target_table for each country and store them separately (one country per csv) into "output_final/202110gavi-3/target_table"
# all countries except COD because it still has some scenarios running
all_countries <- c("AGO", "BDI", "BEN", "BFA", "CAF", "CIV", "CMR", "COD", "COG", "ETH", "GHA", "GIN", "GNB", "KEN", "LBR", "MDG", "MLI", "MOZ", "MRT", "MWI", "NAM", "NER", "NGA", "RWA", "SEN", "SLE", "SOM", "SSD", "TCD", "TGO", "TZA", "UGA", "ZAF", "ZMB", "ZWE")

## make one target table for each country (done for all countries)
for(country in all_countries){
    
    fn <- paste0("target_table_", country, ".csv")
    target_table  <- make_target_table(countries = country, thresholds = c("0.001", "2e-04", "1e-04"),
                                            output_raw_path = "output_raw/202110gavi-3")
    write.csv(target_table, paste0(output_final_path, "/target_table/", fn), row.names = F)
    
}

## make table for targeted population, efficiency and averted cases with CI by country
for(country in all_countries){
    
    # read in target table 
    df_tt <- read.csv(paste0(output_final_path, "/target_table/target_table_", country, ".csv" ))
    df_tt <- df_tt %>% mutate(true_averted_cases = no_vaccination_true_case - campaign_default_true_case)

    # added on 11/10: remove targets smaller than 5*5 grid (averted case == 0 & is_target == 1)
    df_tt <- df_tt %>% filter(true_averted_cases != 0 | is_target != 1)
    
    # get targeted pop and efficiency table for this one country
    df_tp_eff <- combine_eff_table_median(countries = country, df_target = df_tt, admin_level = "both")
    
    # rbind all countries
    if(which(all_countries == country) == 1){df_tp_eff_byISO <- df_tp_eff}
    if(which(all_countries == country) != 1){df_tp_eff_byISO <- rbind(df_tp_eff_byISO, df_tp_eff)}

}

write.csv(df_tp_eff_byISO, paste0(output_final_path, "/eff_table/tp_eff_ac_byISO_median.csv"), row.names = F)


## 10/15 make table for targeted population and efficiency for all countries combined together
# due to memory limit issue, process 5 countries a time
all_countries <- c("AGO", "BDI", "BEN", "BFA", "CAF", 
                   "CIV", "CMR", "COD", "COG", "ETH", 
                   "GHA", "GIN", "GNB", "KEN", "LBR", 
                   "MDG", "MLI", "MOZ", "MRT", "MWI", 
                   "NAM", "NER", "NGA", "RWA", "SEN", 
                   "SLE", "SOM", "SSD", "TCD", "TGO", 
                   "TZA", "UGA", "ZAF", "ZMB", "ZWE")
all_countries <- c("TZA", "UGA", "ZAF", "ZMB", "ZWE")
for(country in all_countries){ # make a big target table stacking all target tables of all countries

    df_tt <- read.csv(paste0(output_final_path, "/target_table/target_table_", country, ".csv" ))
    if(which(all_countries == country) == 1){df_tt_allISOs <- df_tt}
    if(which(all_countries == country) != 1){df_tt_allISOs <- rbind(df_tt_allISOs, df_tt)}

}

df_tt_allISOs <- df_tt_allISOs %>%  mutate(true_averted_cases = no_vaccination_true_case - campaign_default_true_case)

# added on 11/10: remove targets smaller than 5*5 grid (averted case == 0 & is_target == 1)
df_tt_alllISOs <- df_tt_allISOs %>% filter(true_averted_cases != 0 | is_target != 1)

df_tp <- df_tt_allISOs %>%
    #filter(is_target == 1) %>% # added on 10/19 keeping runs even there is no targeting going on
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

df_fvp <- df_tt_allISOs %>%
    group_by(year, run_id, threshold, confirmation_lens, admin_level) %>%
    summarize(fvp = sum(actual_fvp)) %>%
    group_by(run_id, threshold, confirmation_lens, admin_level) %>%
    mutate(fvp_cumu = cumsum(fvp)) %>%
    dplyr::select(-fvp)
    
df_temp <- df_tp %>% 
    full_join(df_ac, by = c("year", "run_id", "threshold", "confirmation_lens", "admin_level")) %>%
    full_join(df_fvp, by = c("year", "run_id", "threshold", "confirmation_lens", "admin_level"))


write.csv(df_temp, paste0(output_final_path, "/temp_table/df_temp_7.csv"), row.names = F)

df_temp <- 
    rbind(
        read.csv(paste0(output_final_path, "/temp_table/df_temp_1.csv")),
        read.csv(paste0(output_final_path, "/temp_table/df_temp_2.csv")),
        read.csv(paste0(output_final_path, "/temp_table/df_temp_3.csv")),
        read.csv(paste0(output_final_path, "/temp_table/df_temp_4.csv")),
        read.csv(paste0(output_final_path, "/temp_table/df_temp_5.csv")),
        read.csv(paste0(output_final_path, "/temp_table/df_temp_6.csv")),
        read.csv(paste0(output_final_path, "/temp_table/df_temp_7.csv"))
    ) %>% 
    group_by(year, run_id, threshold, confirmation_lens, admin_level) %>%
    summarize(target_pop_cumu = sum(target_pop_cumu),
              true_ac_cumu = sum(true_ac_cumu),
              fvp_cumu = sum(fvp_cumu))

write.csv(df_temp, paste0(output_final_path, "/eff_table/tp_ac_fvp_allISOs.csv"), row.names = F)

# calculate median and CI (all ISOs combined)
df_tp <- df_temp %>%
    group_by(year, threshold, confirmation_lens, admin_level) %>%
    #summarize(tp_cumu_mean = mean(target_pop_cumu, na.rm = TRUE),
    #          tp_cumu_sd = sd(target_pop_cumu, na.rm = TRUE),
    #          tp_cumu_n = n()) %>%
    #mutate(tp_cumu_se = tp_cumu_sd / sqrt(tp_cumu_n),
    #       tp_cumu_lb = tp_cumu_mean - qt(1 - (0.05 / 2), tp_cumu_n - 1) * tp_cumu_se,
    #       tp_cumu_ub = tp_cumu_mean + qt(1 - (0.05 / 2), tp_cumu_n - 1) * tp_cumu_se) %>%
    summarize(tp_cumu_median = median(target_pop_cumu, na.rm = T),
              tp_cumu_lb = quantile(target_pop_cumu, 0.025, na.rm = T),
              tp_cumu_ub = quantile(target_pop_cumu, 0.975, na.rm = T)) %>%
    dplyr::select(year, threshold, confirmation_lens, admin_level, tp_cumu_lb, tp_cumu_median, tp_cumu_ub) %>%
    group_by(threshold, confirmation_lens, admin_level) %>%
    summarize(tp_cumu_lb = max(tp_cumu_lb, na.rm = T),
              tp_cumu_median = max(tp_cumu_median, na.rm = T),
              tp_cumu_ub = max(tp_cumu_ub, na.rm = T)) %>%
    dplyr::select(threshold, confirmation_lens, admin_level, tp_cumu_lb, tp_cumu_median, tp_cumu_ub)
    
df_eff <- df_temp %>%
   filter(year == 2035) %>%
   mutate(efficiency = true_ac_cumu / fvp_cumu * 1000) %>%
   group_by(threshold, confirmation_lens, admin_level) %>%
   #summarize(efficiency_mean = mean(efficiency, na.rm = TRUE),
   #          efficiency_sd = sd(efficiency, na.rm = TRUE),
   #          efficiency_n = n()) %>%
   #mutate(efficiency_se = efficiency_sd / sqrt(efficiency_n),
   #       efficiency_lb = efficiency_mean - qt(1 - (0.05 / 2), efficiency_n - 1) * efficiency_se,
   #       efficiency_ub = efficiency_mean + qt(1 - (0.05 / 2), efficiency_n - 1) * efficiency_se) %>%
   summarize(efficiency_median = median(efficiency, na.rm = TRUE),
             efficiency_lb = quantile(efficiency, 0.025, na.rm = TRUE),
             efficiency_ub = quantile(efficiency, 0.975, na.rm = TRUE)) %>%
   dplyr::select(threshold, confirmation_lens, admin_level, efficiency_lb, efficiency_median, efficiency_ub)


df_ac <- df_temp %>%
    filter(year == 2035) %>%
    group_by(threshold, confirmation_lens, admin_level) %>%
    #summarize(ac_mean = mean(true_ac_cumu, na.rm = TRUE),
    #            ac_sd = sd(true_ac_cumu, na.rm = TRUE),
    #            ac_n = n()) %>%
    #mutate(ac_se = ac_sd / sqrt(ac_n),
    #       ac_lb = ac_mean - qt(1 - (0.05 / 2), ac_n - 1) * ac_se,
    #       ac_ub = ac_mean + qt(1 - (0.05 / 2), ac_n - 1) * ac_se) %>%
    summarize(ac_median = median(true_ac_cumu, na.rm = TRUE),
              ac_lb = quantile(true_ac_cumu, 0.025, na.rm = TRUE),
              ac_ub = quantile(true_ac_cumu, 0.975, na.rm = TRUE)) %>%
    dplyr::select(threshold, confirmation_lens, admin_level, ac_lb, ac_median, ac_ub)


# join two tables and round up numbers
df_result <- df_tp %>%
    full_join(df_eff, by = c("threshold", "confirmation_lens", "admin_level")) %>%
    full_join(df_ac, by = c("threshold", "confirmation_lens", "admin_level"))
    # convert unit of tp to "million"
    # mutate(tp_cumu_lb = round(tp_cumu_lb/1000000, 2),
    #         tp_cumu_median = round(tp_cumu_median/1000000, 2),
    #         tp_cumu_ub = round(tp_cumu_ub/1000000, 2),
    #         efficiency_lb = round(efficiency_lb, 2),
    #         efficiency_median = round(efficiency_median, 2),
    #         efficiency_ub = round(efficiency_ub, 2)) %>%
    #  mutate(tp_cumu = paste0(tp_cumu_median, " (", tp_cumu_lb, "-", tp_cumu_ub, ")"),
    #         efficiency = paste0(efficiency_median, " (", efficiency_lb, "-", efficiency_ub, ")")) %>%
    #  select(threshold, confirmation_lens, admin_level, tp_cumu, efficiency)

write.csv(df_result, paste0(output_final_path, "/eff_table/tp_eff_ac_allISOs_median.csv"), row.names = F)


## added on 11/9 to get number of targets of all years
# aim to get numbers (across all years) of 1) total targets and 2) unique targets:
library(dplyr)
library(stringr)
output_final_path <- "output_final/202110gavi-3"
all_countries <- c("AGO", "BDI", "BEN", "BFA", "CAF", "CIV", "CMR", "COD", "COG", "ETH", "GHA", "GIN", "GNB", "KEN", "LBR", "MDG", "MLI", "MOZ", "MRT", "MWI", "NAM", "NER", "NGA", "RWA", "SEN", "SLE", "SOM", "SSD", "TCD", "TGO", "TZA", "UGA", "ZAF", "ZMB", "ZWE")

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

n_targeted_admins <- 
    full_join(n_total_targets, n_unique_targets,
              by = c("ISO", "confirmation_lens", "admin_level", "threshold"))

write.csv(n_targeted_admins, "output_final/202110gavi-3/eff_table/n_targeted_admins_byISO.csv", row.names = F)

# combine all countries
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

n_targeted_admins <- 
    full_join(n_total_targets, n_unique_targets,
              by = c("confirmation_lens", "admin_level", "threshold"))

write.csv(n_targeted_admins, "output_final/202110gavi-3/eff_table/n_targeted_admins_allISOs.csv", row.names = F)


# added on 11/9/22 --------------------------------------------------
# to look at how many of places are missed by the fasterize() function
library(dplyr)
library(stringr)
output_final_path <- "output_final/202110gavi-3"
all_countries <- c("AGO", "BDI", "BEN", "BFA", "CAF", "CIV", "CMR", "COD", "COG", "ETH", "GHA", "GIN", "GNB", "KEN", "LBR", "MDG", "MLI", "MOZ", "MRT", "MWI", "NAM", "NER", "NGA", "RWA", "SEN", "SLE", "SOM", "SSD", "TCD", "TGO", "TZA", "UGA", "ZAF", "ZMB", "ZWE")

for(country in all_countries){
    
    df_tt <- read.csv(paste0(output_final_path, "/target_table/target_table_", country, ".csv" ))
    df_tt <- df_tt %>% mutate(true_averted_cases = no_vaccination_true_case - campaign_default_true_case)
    temp <- filter(df_tt, is_target == 1 & true_averted_cases == 0)
    
    if(which(all_countries == country) == 1){df_problem <- temp}
    if(which(all_countries == country) != 1){df_problem <- rbind(df_problem, temp)}

}

write.csv(df_problem, "output_final/202110gavi-3/debug/df_small_targets.csv", row.names = F)


# added on 11/10: numbers runs with targeting going on
for(country in all_countries){
    
    df_tt <- read.csv(paste0(output_final_path, "/target_table/target_table_", country, ".csv" ))
    df_tt <- df_tt %>% mutate(true_averted_cases = no_vaccination_true_case - campaign_default_true_case)
    # remove those small districts that were targeted
    df_tt <- df_tt %>% filter(is_target != 1 | true_averted_cases != 0)
    
    n_targeted_runs <- df_tt %>%
        group_by(ISO, confirmation_lens, threshold, admin_level, run_id) %>%
        summarize(sum_target = sum(is_target)) %>%
        group_by(ISO, confirmation_lens, threshold, admin_level) %>%
        count(sum_target != 0) %>%
        filter(`sum_target != 0` == TRUE) %>%
        rename(n_runs_targeted = n) %>%
        dplyr::select(ISO, threshold, confirmation_lens, admin_level, n_runs_targeted)

    if(which(all_countries == country) == 1){df_n_targeted_runs <- n_targeted_runs}
    if(which(all_countries == country) != 1){df_n_targeted_runs <- rbind(df_n_targeted_runs, n_targeted_runs)}
    
}

write.csv(df_n_targeted_runs, "output_final/202110gavi-3/eff_table/n_targeted_runs.csv", row.names = F)



## (11/16) Look into what might be the main reason for high variability in 10/10,000 scenario

# select one country as an example 
country <- "SOM"
df_tt <- read.csv(paste0(output_final_path, "/target_table/target_table_", country, ".csv" ))
df_tt <- df_tt %>% mutate(true_averted_cases = no_vaccination_true_case - campaign_default_true_case)
df_tt <- df_tt %>% filter(threshold == 1e-04) %>% filter(is_target == 1)

# Two possible scenarios to cause the variability:
# 1. the incidence rate raster are on different ends, when incidence rates happen to be at lower end - targets a lot of places, vice verse;
# 2. has many variability in the incidence rate  - different targets every layer

# look at number of targets, and targets of in a certain year, under certain threshold 
temp <- df_tt %>% filter(year == 2022) %>%
    filter(confirmation_lens == "district-estimate" & admin_level == "admin2") %>%
    group_by(run_id) %>%
    arrange(confirmed_incidence_rate) %>%
    summarize(n = n(), targets = paste(NAME_2, collapse = ", "))

write.csv(temp, paste0("targets_2022_", country, ".csv"), row.names = F)


## 11/17 Answering Elizabeth's comment:
# Could I see "proportion of country-simulations with no targets" and 
# associated number of unique countries with no targets by scenario? 
# Essentially collapsing the country stratification for that Table S2?
setwd("/home/hanmeng/VIMC/gavi_vimc_cholera")
library(tidyr)
all_countries <- c("AGO", "BDI", "BEN", "BFA", "CAF", "CIV", "CMR", "COD", "COG", "ETH", "GHA", "GIN", "GNB", "KEN", "LBR", "MDG", "MLI", "MOZ", "MRT", "MWI", "NAM", "NER", "NGA", "RWA", "SEN", "SLE", "SOM", "SSD", "TCD", "TGO", "TZA", "UGA", "ZAF", "ZMB", "ZWE")

df_result <- expand_grid(
    threshold = c(1e-04, 2e-04, 1e-03),
    confirmation_lens = c("district-estimate", "global-estimate", "no-estimate"),
    admin_level = c("admin1", "admin2")
)

## look at whether a country has any targeting going on by scenario
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
## if there is no targeting going on in this country for all the simulations (for a given scenario), then noted as 0, otherwise 1.
write.csv(df_targeted, paste0(output_final_path, "/eff_table/ind_targeted.csv"), row.names = FALSE)


## Then look at for each simulation, whether a country has been targeted
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
        select(-n)
    colnames(temp)[5] <- country
    
    if(which(all_countries == country) == 1){df_targeted_byrun <- temp}
    if(which(all_countries == country) != 1){df_targeted_byrun <- df_targeted_byrun %>% left_join(temp,  by = c("run_id", "threshold", "confirmation_lens", "admin_level"))}

}
write.csv(df_targeted_byrun, paste0(output_final_path, "/eff_table/ind_targeted_byrun.csv"), row.names = FALSE)


## (11/22) Number of unique countries targeted
setwd(paste0(output_final_path, "/eff_table"))
ind_targeted_byrun <- read.csv("ind_targeted_byrun.csv")
n_targeted_ISOs <- ind_targeted_byrun %>%
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
n_targeted_admins_allISOs <- read.csv("n_targeted_admins_allISOs.csv")
n_targeted_admins_ISOs <- n_targeted_admins_allISOs %>%
    left_join(n_targeted_ISOs, by = c("confirmation_lens", "threshold", "admin_level"))

write.csv(n_targeted_admins_ISOs, "n_targeted_admins_ISOs.csv", row.names = FALSE)