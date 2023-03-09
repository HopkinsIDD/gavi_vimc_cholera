# Script to summarize cost and testing metrics

# set up =================================
r_lib <- Sys.getenv("R_LIBRARY_DIRECTORY", FALSE)
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
source("packages/ocvImpact/R/utils_cost.R")
output_final_path <- "output_final/202110gavi-3"
all_countries <- c("AGO", "BDI", "BEN", "BFA", "CAF", "CIV", "CMR", "COD", "COG", "ETH", "GHA", "GIN", "GNB", "KEN", "LBR", "MDG", "MLI", "MOZ", "MRT", "MWI", "NAM", "NER", "NGA", "RWA", "SEN", "SLE", "SOM", "SSD", "TCD", "TGO", "TZA", "UGA", "ZAF", "ZMB", "ZWE")


# set parameters (exact values pending) =================
cost_per_ocv_procure <- 1.65
cost_per_ocv_deliver <- 0.65
cost_per_ocv_shipping <- 0.1
cost_per_ocv <- cost_per_ocv_procure + cost_per_ocv_deliver + cost_per_ocv_shipping
cost_per_rdt <- 1.9
cost_per_culture <- 12.84
rdt_pos_rate <- 0.56
prop_sus_rdt_tested <- 0.74 # proportion of suspected (clinical) cases tested by RDT at district confirmation capacity
prop_rdt_positive_cultured <- 0.48 # proportion of RDT+ cases tested by culture at district confirmation capacity
prop_sus_cultured <- 0.3 # proportion of suspected (clinical cases) tested by culture at national confirmation capacity
dose_regimen <- 2


# check directory status =====================
if(!file.exists(paste0(output_final_path, "/", "cost_table"))){
    dir.create(paste0(output_final_path, "/", "cost_table"))
}else{
    if(length(dir(all.files=TRUE)) != 0){
        message("Directory output_final/202110gavi-3/cost_table is not empty and cost tables might have already been made.")
    }
}

# set output paths ===========
int_path <- paste0(output_final_path, "/intermediate_table/")
cost_path <- paste0(output_final_path, "/cost_table/")
tt_path <- paste0(output_final_path, "/target_table/")


# ===============================================
# PREP: Extract clinical cases from target tables 
# ===============================================
# OUTPUT 13A clinical cases by country
for(country in all_countries){
    
    # read in target table 
    message(paste0("For making dataframe of clinical cases by country, reading in target table of ", country))
    df_tt <- read.csv(paste0(tt_path, "target_table_", country, ".csv" ))
    df_tt <- df_tt %>% mutate(true_averted_cases = no_vaccination_true_case - campaign_default_true_case)

    # added on 11/10/2022: remove targets smaller than 5*5 grid (averted case == 0 & is_target == 1)
    df_tt <- df_tt %>% filter(true_averted_cases != 0 | is_target != 1)
    
    # OUTPUT 13A: clinical cases by country with median and CI (district-estimate and global-estimate only)
    # get clinical cases for district-estimate and global-estimate scenarios (assuming no testing in no-vaccination and no-estimate scenario)
    df_clinical_medianCI <- get_clinical_cases_medianCI(df_tt)
    message(paste0("Rbinding the clinical cases of ", country, " to existing data clinical_cases_byISO_medianCI data frame."))
    if(which(all_countries == country) == 1){clinical_cases_byISO_medianCI <- df_clinical_medianCI}
    if(which(all_countries == country) != 1){clinical_cases_byISO_medianCI <- rbind(clinical_cases_byISO_medianCI, df_clinical_medianCI)}
   
}
rm(df_tt, df_clinical_medianCI)
write.csv(clinical_cases_byISO_medianCI, paste0(int_path, "clinical_cases_byISO_medianCI.csv"), row.names = F)


# OUTPUT 13B clinical cases across countries
for(country in all_countries){
    
    # read in target table 
    message(paste0("For making dataframe of clinical cases by country and by run, reading in target table of ", country))
    df_tt <- read.csv(paste0(tt_path, "target_table_", country, ".csv" ))
    df_tt <- df_tt %>% mutate(true_averted_cases = no_vaccination_true_case - campaign_default_true_case)

    # added on 11/10/2022: remove targets smaller than 5*5 grid (averted case == 0 & is_target == 1)
    df_tt <- df_tt %>% filter(true_averted_cases != 0 | is_target != 1)
    
    # OUTPUT 13B: clinical cases across all countries over the years by run (large memory requirement)
    df_clinical <- get_clinical_cases(df_tt)
    message(paste0("Rbinding the clinical cases of ", country, " to existing data df_clinical_byISO data frame."))
    if(which(all_countries == country) == 1){df_clinical_byISO <- df_clinical}
    if(which(all_countries == country) != 1){df_clinical_byISO <- rbind(df_clinical_byISO, df_clinical)}

}
rm(df_tt, df_clinical)
clinical_cases_allISOs <- df_clinical_byISO %>%
    group_by(confirmation_lens, admin_level, threshold, run_id) %>%
    summarize(cumu_clinical_cases = sum(cumu_clinical_cases))
write.csv(clinical_cases_allISOs, paste0(int_path, "clinical_cases_allISOs.csv"), row.names = F)



# =====================================
# PART 1: Cost of vaccination campaigns 
# =====================================
# OUTPUT 14: ocv_cost_byISO_medianCI -----------------------
# OCV campaign cost for each country, with median and CI
if(!file.exists(paste0(int_path, "tp_eff_ac_byISO_medianCI.csv"))){
    message("tp_eff_ac_byISO_medianCI.csv doesn't exist in the intermediate_table folder, please check." )
}else{
    tp_eff_ac_byISO_medianCI <- read.csv(paste0(int_path, "tp_eff_ac_byISO_medianCI.csv"))
}

ocv_cost_byISO_medianCI <- tp_eff_ac_byISO_medianCI %>%
    mutate(ocv_cost_lb = dose_regimen * cost_per_ocv * tp_cumu_lb,
           ocv_cost_median = dose_regimen * cost_per_ocv * tp_cumu_median,
           ocv_cost_ub = dose_regimen * cost_per_ocv * tp_cumu_ub) %>%
    dplyr::select(ISO, threshold, confirmation_lens, admin_level, 
                  ocv_cost_lb, ocv_cost_median, ocv_cost_ub)

write.csv(ocv_cost_byISO_medianCI, paste0(cost_path, "ocv_cost_byISO_medianCI.csv"), row.names = F)

# OUTPUT 15: ocv_cost_allISOs_medianCI -----------------------
# OCV campaign cost all countries combined, with median and CI
if(!file.exists(paste0(int_path, "tp_eff_ac_allISOs_medianCI.csv"))){
    message("tp_eff_ac_allISOs_medianCI.csv doesn't exist in the intermediate_table folder, please check." )
}else{
    tp_eff_ac_allISOs_medianCI <- read.csv(paste0(int_path, "tp_eff_ac_allISOs_medianCI.csv"))
}

ocv_cost_allISOs_medianCI <- tp_eff_ac_allISOs_medianCI %>%
    mutate(ocv_cost_lb = tp_cumu_lb * dose_regimen * cost_per_ocv,
           ocv_cost_median = tp_cumu_median * dose_regimen * cost_per_ocv,
           ocv_cost_ub = tp_cumu_ub * dose_regimen * cost_per_ocv) %>%
    dplyr::select(threshold, confirmation_lens, admin_level, 
                  ocv_cost_lb, ocv_cost_median, ocv_cost_ub)

write.csv(ocv_cost_allISOs_medianCI, paste0(cost_path, "ocv_cost_allISOs_medianCI.csv"), row.names = F)

# =====================================
# PART 2: Cost of testing (culture, RDT)
# =====================================
# OUTPUT 16: test_cost_byISO_medianCI ---------------------
if(!file.exists(paste0(int_path, "clinical_cases_byISO_medianCI.csv"))){
    message("Table of clinical cases (w/ median and CI) doesn't exist in the intermediate_table folder, please check." )
}else{
    clinical_cases_byISO_medianCI <- read.csv(paste0(int_path, "clinical_cases_byISO_medianCI.csv"))
}

test_cost_byISO_medianCI <- clinical_cases_byISO_medianCI %>%
    # RDT test cost
    mutate(rdt_cost_lb = case_when(confirmation_lens == "global-estimate" ~ 0,
                                   confirmation_lens == "district-estimate" ~ cumu_clinical_cases_lb * prop_sus_rdt_tested * cost_per_rdt)) %>%
    mutate(rdt_cost_median = case_when(confirmation_lens == "global-estimate" ~ 0,
                                       confirmation_lens == "district-estimate" ~ cumu_clinical_cases_median * prop_sus_rdt_tested * cost_per_rdt)) %>%
    mutate(rdt_cost_ub = case_when(confirmation_lens == "global-estimate" ~ 0,
                                   confirmation_lens == "district-estimate" ~ cumu_clinical_cases_ub * prop_sus_rdt_tested * cost_per_rdt)) %>%
    # culture test cost 
    mutate(culture_cost_lb = case_when(confirmation_lens == "global-estimate" ~ cumu_clinical_cases_lb * prop_sus_cultured * cost_per_culture,
                                       confirmation_lens == "district-estimate" ~ cumu_clinical_cases_lb * rdt_pos_rate * prop_rdt_positive_cultured * cost_per_culture)) %>%
    mutate(culture_cost_median = case_when(confirmation_lens == "global-estimate" ~ cumu_clinical_cases_median * prop_sus_cultured * cost_per_culture,
                                           confirmation_lens == "district-estimate" ~ cumu_clinical_cases_median * rdt_pos_rate * prop_rdt_positive_cultured * cost_per_culture)) %>%
    mutate(culture_cost_ub = case_when(confirmation_lens == "global-estimate" ~ cumu_clinical_cases_ub * prop_sus_cultured * cost_per_culture,
                                       confirmation_lens == "district-estimate" ~ cumu_clinical_cases_ub * rdt_pos_rate * prop_rdt_positive_cultured * cost_per_culture)) %>%
    # Total (RDT + culture) test cost 
    mutate(test_cost_lb = rdt_cost_lb + culture_cost_lb,
           test_cost_median = rdt_cost_median + culture_cost_median,
           test_cost_ub = rdt_cost_ub + culture_cost_ub) %>%
    dplyr::select(-cumu_clinical_cases_lb, -cumu_clinical_cases_median, -cumu_clinical_cases_ub)

write.csv(test_cost_byISO_medianCI, paste0(cost_path, "test_cost_byISO_medianCI.csv"), row.names = F)


# OUTPUT 17: test_cost_allISOs_medianCI ---------------------
if(!file.exists(paste0(cost_path, "test_cost_byISO_medianCI.csv"))){
    message("Dataframe test_cost_byISO_medianCI.csv doesn't exist in the intermediate_table folder, please check." )
}else{
    test_cost_byISO_medianCI <- read.csv(paste0(cost_path, "test_cost_byISO_medianCI.csv"))
}

test_cost_allISOs_medianCI <- test_cost_byISO_medianCI %>%
    group_by(confirmation_lens, admin_level, threshold) %>%
    summarize(rdt_cost_lb = sum(rdt_cost_lb),
              rdt_cost_median = sum(rdt_cost_median),
              rdt_cost_ub = sum(rdt_cost_ub),
              culture_cost_lb = sum(culture_cost_lb),
              culture_cost_median = sum(culture_cost_median),
              culture_cost_ub = sum(culture_cost_ub),
              test_cost_lb = sum(test_cost_lb),
              test_cost_median = sum(test_cost_median),
              test_cost_ub = sum(test_cost_ub))

write.csv(test_cost_allISOs_medianCI, paste0(cost_path, "test_cost_allISOs_medianCI.csv"), row.names = F)


# ===================================================
# PART 3: Cost-effectiveness (vaccine campaigns only) 
# ===================================================
# OUTPUT 18: ocv_costeff_byISO_medianCI 
# cost effectiveness of OCV campaign by country
ocv_costeff_byISO_medianCI <- tp_eff_ac_byISO_medianCI %>%
    mutate(ocv_costeff_lb = 1000 * dose_regimen * cost_per_ocv / efficiency_ub,
           ocv_costeff_median = 1000 * dose_regimen * cost_per_ocv / efficiency_median,
           ocv_costeff_ub = 1000 * dose_regimen * cost_per_ocv / efficiency_lb) %>%
    dplyr::select(ISO, threshold, confirmation_lens, admin_level, 
                  ocv_costeff_lb, ocv_costeff_median, ocv_costeff_ub)

write.csv(ocv_costeff_byISO_medianCI, paste0(cost_path, "ocv_costeff_byISO_medianCI.csv"), row.names = F)

# OUTPUT 19: ocv_costeff_allISOs_medianCI 
ocv_costeff_allISOs_medianCI <- tp_eff_ac_allISOs_medianCI %>%
    mutate(ocv_costeff_lb = 1000 * dose_regimen * cost_per_ocv / efficiency_ub,
           ocv_costeff_median = 1000 * dose_regimen * cost_per_ocv / efficiency_median,
           ocv_costeff_ub = 1000 * dose_regimen * cost_per_ocv / efficiency_lb) %>%
    dplyr::select(threshold, confirmation_lens, admin_level, 
                  ocv_costeff_lb, ocv_costeff_median, ocv_costeff_ub)

write.csv(ocv_costeff_allISOs_medianCI, paste0(cost_path, "ocv_costeff_allISOs_medianCI.csv"), row.names = F)

# ========================================================
# PART 4: Cost-effectiveness (vaccine campaigns + testing) 
# ========================================================
# OUTPUT 20: total_costeff_byISO_medianCI 
# Total cost effectiveness (vaccine campaigns + testing) 
# SEE ABOVE (already combined in the same for loop to reduce run time)
for(country in all_countries){
    
    # read in target table 
    message(paste0("For making total_costeff_byISO_medianCI, reading target table of ", country))
    df_tt <- read.csv(paste0(tt_path, "target_table_", country, ".csv" ))
    df_tt <- df_tt %>% mutate(true_averted_cases = no_vaccination_true_case - campaign_default_true_case)

    # added on 11/10/2022: remove targets smaller than 5*5 grid (averted case == 0 & is_target == 1)
    df_tt <- df_tt %>% filter(true_averted_cases != 0 | is_target != 1)
    
    # get clinical cases for district-estimate and global-estimate scenarios (assuming no testing in no-vaccination and no-estimate scenario)
    df_cost <- calc_costeff_ocv_test(df_tt)
    
    rm(df_tt)
    message(paste0("Rbinding costeff of ", country, " to existing total_costeff_byISO_medianCI."))
    if(which(all_countries == country) == 1){total_costeff_byISO_medianCI <- df_cost}
    if(which(all_countries == country) != 1){total_costeff_byISO_medianCI <- rbind(total_costeff_byISO_medianCI, df_cost)}
    rm(df_cost)
}

write.csv(total_costeff_byISO_medianCI, paste0(cost_path, "total_costeff_byISO_medianCI.csv"), row.names = F)


# OUTPUT 21: total_costeff_allISOs_medianCI 
# total cost effectiveness (ocv campaign + test)
if(!file.exists(paste0(int_path, "tp_ac_allISOs.csv")) |
   !file.exists(paste0(int_path, "clinical_cases_allISOs.csv"))){
    message("At least one of tp_ac_allISOs.csv or clinical_cases_allISOs.csv doesn't exist in the intermediate_table folder, please check." )
}else{
    tp_ac_allISOs <- read.csv(paste0(int_path, "tp_ac_allISOs.csv"))
    clinical_cases_allISOs <- read.csv(paste0(int_path, "clinical_cases_allISOs.csv"))
}


total_costeff_allISOs_medianCI <- tp_ac_allISOs %>%
    filter(year == 2035) %>% dplyr::select(-year) %>%
    left_join(clinical_cases_allISOs, by = c("confirmation_lens", "admin_level", "threshold", "run_id")) %>%
    mutate(ocv_cost = target_pop_cumu * dose_regimen * cost_per_ocv) %>%
    mutate(rdt_cost = case_when(confirmation_lens == "global-estimate" | confirmation_lens == "no-estimate" ~ 0,
                                confirmation_lens == "district-estimate" ~ cumu_clinical_cases * prop_sus_rdt_tested * cost_per_rdt)) %>%
    mutate(culture_cost = case_when(confirmation_lens == "no-estimate" ~ 0,
                                    confirmation_lens == "global-estimate" ~ cumu_clinical_cases * prop_sus_cultured * cost_per_culture,
                                    confirmation_lens == "district-estimate" ~ cumu_clinical_cases * rdt_pos_rate * prop_rdt_positive_cultured * cost_per_culture)) %>%
    mutate(test_cost = rdt_cost + culture_cost) %>%
    mutate(costeff_ocv_test = (ocv_cost + test_cost) / true_ac_cumu) %>% 
    group_by(confirmation_lens, admin_level, threshold) %>%
    summarize(costeff_ocv_test_lb = quantile(costeff_ocv_test, 0.025, na.rm = T),
              costeff_ocv_test_median = quantile(costeff_ocv_test, 0.5, na.rm = T),
              costeff_ocv_test_ub = quantile(costeff_ocv_test, 0.975, na.rm = T)) 
        
write.csv(total_costeff_allISOs_medianCI, paste0(cost_path, "total_costeff_allISOs_medianCI.csv"), row.names = F)


# ==============================
# PART 5: ICER between scenarios
# ==============================
# ICER (Incremental Cost-Effectiveness Ratio): the ratio of the difference in costs between two strategies to the difference in effectiveness (averted cases).
# Two sets of comparison:
    # 1) compare with no vaccination scenario: no cost at all.
	#    ICER = (cost of scenario X)/(averted case of scenario X)) 
	#    ----> this is the same as total_costeff calculated above.
    # 2) compare with "status quo" scenario: only vaccination cost
	#    ICER = (cost of scenario X  - cost of status quo scenario) /
	#           (averted cases of scenario X - averted cases of status quo scenario)
    # (status quo scenario: clinical case detection + district targeting + 2 per 10,000 incidence rate threshold.)
status_quo_confirmation_lens <- "no-estimate"
status_quo_admin_level <- "admin2"
status_quo_threshold <- 2e-04 # should be 2e-04

# OUTPUT 22: ICER_allISOs (ICER summarizing all countries)
if(!file.exists(paste0(int_path, "tp_ac_allISOs.csv")) |
   !file.exists(paste0(int_path, "clinical_cases_allISOs.csv"))){
    message("At least one of tp_ac_allISOs.csv or clinical_cases_allISOs.csv doesn't exist in the intermediate_table folder, please check." )
}else{
    tp_ac_allISOs <- read.csv(paste0(int_path, "tp_ac_allISOs.csv"))
    clinical_cases_allISOs <- read.csv(paste0(int_path, "clinical_cases_allISOs.csv"))
}

df_cost_ac <- tp_ac_allISOs %>%
    filter(year == 2035) %>% dplyr::select(-year) %>%
    left_join(clinical_cases_allISOs, by = c("confirmation_lens", "admin_level", "threshold", "run_id")) %>%
    mutate(ocv_cost = target_pop_cumu * dose_regimen * cost_per_ocv) %>%
    mutate(rdt_cost = case_when(confirmation_lens == "global-estimate" | confirmation_lens == "no-estimate" ~ 0,
                                confirmation_lens == "district-estimate" ~ cumu_clinical_cases * prop_sus_rdt_tested * cost_per_rdt)) %>%
    mutate(culture_cost = case_when(confirmation_lens == "no-estimate" ~ 0,
                                    confirmation_lens == "global-estimate" ~ cumu_clinical_cases * prop_sus_cultured * cost_per_culture,
                                    confirmation_lens == "district-estimate" ~ cumu_clinical_cases * rdt_pos_rate * prop_rdt_positive_cultured * cost_per_culture)) %>%
    mutate(test_cost = rdt_cost + culture_cost) %>%
    mutate(total_cost = ocv_cost + test_cost)

df_icer <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(df_icer) <- c("run_id", "threshold", "confirmation_lens", "admin_level", "icer")

for(confirmation_lens_selected in c("no-estimate", "global-estimate", "district-estimate")){
    for(admin_level_selected in c("admin1", "admin2")){
        for(threshold_selected in unique(df_cost_ac$threshold)){
            df_icer_temp <- df_cost_ac %>%
                dplyr::select(run_id, threshold, confirmation_lens, admin_level, total_cost, true_ac_cumu) %>%
                group_by(run_id) %>%
                mutate(icer = (total_cost[confirmation_lens == confirmation_lens_selected & admin_level == admin_level_selected & threshold == threshold_selected] - total_cost[confirmation_lens == status_quo_confirmation_lens & admin_level == status_quo_admin_level & threshold == status_quo_threshold]) /
                              (true_ac_cumu[confirmation_lens == confirmation_lens_selected & admin_level == admin_level_selected & threshold == threshold_selected] - true_ac_cumu[confirmation_lens == status_quo_confirmation_lens & admin_level == status_quo_admin_level & threshold == status_quo_threshold])) %>%
                filter(confirmation_lens == confirmation_lens_selected & admin_level == admin_level_selected & threshold == threshold_selected) %>%
                dplyr::select(-total_cost, -true_ac_cumu)
                
            df_icer <- rbind(df_icer, df_icer_temp)
        }
    }
}

# get median and CI of ICER
icer_allISOs_medianCI <- df_icer %>%
    group_by(threshold, confirmation_lens, admin_level) %>%
    summarize(icer_lb = quantile(icer, 0.025, na.rm = T),
              icer_median = quantile(icer, 0.5, na.rm = T),
              icer_ub = quantile(icer, 0.975, na.rm = T))

write.csv(icer_allISOs_medianCI, paste0(cost_path, "icer_allISOs_medianCI.csv"), row.names = F)


# OUTPUT 23: icer_byISO_medianCI (ICER of each country)
# get the ICER (with median and CI) of each country
for(country in all_countries){
    
    # read in target table 
    message(paste0("For ICER calculation, reading in target table of ", country))
    df_tt <- read.csv(paste0(tt_path, "target_table_", country, ".csv" ))
    df_tt <- df_tt %>% mutate(true_averted_cases = no_vaccination_true_case - campaign_default_true_case)

    # added on 11/10/2022: remove targets smaller than 5*5 grid (averted case == 0 & is_target == 1)
    df_tt <- df_tt %>% filter(true_averted_cases != 0 | is_target != 1)

    df_tp_ac_temp <- df_tt %>%
        group_by(run_id, year, confirmation_lens, admin_level, threshold) %>%
        summarize(target_pop = sum(actual_fvp), ac = sum(true_averted_cases)) %>%
        group_by(run_id, confirmation_lens, admin_level, threshold) %>%
        mutate(target_pop_cumu = cumsum(target_pop), 
               ac_cumu = cumsum(ac)) %>%
        filter(year == 2035) %>% dplyr::select(-c(year, target_pop, ac))
    
    df_clinical_temp <- get_clinical_cases(df_tt)

    df_tp_ac_clinical_temp <- df_clinical_temp %>%
        left_join(df_tp_ac_temp, by = c("run_id", "confirmation_lens", "threshold", "admin_level"))
    
    # calculate ICER for this individual country
    df_cost_ac_temp <- df_tp_ac_clinical_temp %>%
        mutate(ocv_cost = target_pop_cumu * dose_regimen * cost_per_ocv) %>%
        mutate(rdt_cost = case_when(confirmation_lens == "global-estimate" | confirmation_lens == "no-estimate" ~ 0,
                                    confirmation_lens == "district-estimate" ~ cumu_clinical_cases * prop_sus_rdt_tested * cost_per_rdt)) %>%
        mutate(culture_cost = case_when(confirmation_lens == "no-estimate" ~ 0,
                                        confirmation_lens == "global-estimate" ~ cumu_clinical_cases * prop_sus_cultured * cost_per_culture,
                                        confirmation_lens == "district-estimate" ~ cumu_clinical_cases * rdt_pos_rate * prop_rdt_positive_cultured * cost_per_culture)) %>%
        mutate(test_cost = rdt_cost + culture_cost) %>%
        mutate(total_cost = ocv_cost + test_cost)
    
    df_icer <- data.frame(matrix(ncol = 5, nrow = 0))
    colnames(df_icer) <- c("run_id", "threshold", "confirmation_lens", "admin_level", "icer")

    for(confirmation_lens_selected in c("no-estimate", "global-estimate", "district-estimate")){
        for(admin_level_selected in c("admin1", "admin2")){
            for(threshold_selected in unique(df_cost_ac_temp$threshold)){
                df_icer_temp <- df_cost_ac_temp %>%
                    dplyr::select(run_id, threshold, confirmation_lens, admin_level, total_cost, ac_cumu) %>%
                    group_by(run_id) %>%
                    mutate(icer = (total_cost[confirmation_lens == confirmation_lens_selected & admin_level == admin_level_selected & threshold == threshold_selected] - total_cost[confirmation_lens == status_quo_confirmation_lens & admin_level == status_quo_admin_level & threshold == status_quo_threshold]) /
                                  (ac_cumu[confirmation_lens == confirmation_lens_selected & admin_level == admin_level_selected & threshold == threshold_selected] - ac_cumu[confirmation_lens == status_quo_confirmation_lens & admin_level == status_quo_admin_level & threshold == status_quo_threshold])) %>%
                    filter(confirmation_lens == confirmation_lens_selected & admin_level == admin_level_selected & threshold == threshold_selected) %>%
                    dplyr::select(-total_cost, -ac_cumu)
                
                df_icer <- rbind(df_icer, df_icer_temp)
            }
        }
    }
    
    # get median and CI of ICER of this individual country  
    icer_byISO_medianCI_temp <- df_icer %>%
        group_by(threshold, confirmation_lens, admin_level) %>%
        summarize(icer_lb = quantile(icer, 0.025, na.rm = T),
                  icer_median = quantile(icer, 0.5, na.rm = T),
                  icer_ub = quantile(icer, 0.975, na.rm = T)) %>%
        mutate(ISO = country)
        
    rm(df_tt, df_tp_ac_temp, df_clinical_temp,  df_tp_ac_clinical_temp, df_cost_ac_temp, df_icer_temp, df_icer)
    
    message(paste0("Rbinding ICER table of ", country, " to existing icer_byISO_medianCI."))
    if(which(all_countries == country) == 1){icer_byISO_medianCI <- icer_byISO_medianCI_temp }
    if(which(all_countries == country) != 1){icer_byISO_medianCI <- rbind(icer_byISO_medianCI, icer_byISO_medianCI_temp)}
}
rm(icer_byISO_medianCI_temp)
write(icer_byISO_medianCI, paste0(cost_path, "icer_byISO_medianCI.csv"), row.names = F)


## OUTPUT 28: n_tested_byISO 
# Number of clinical cases that were tested by RDT and/or by culture
if(!file.exists(paste0(output_final_path, "/intermediate_table/clinical_cases_byISO_medianCI.csv"))){
    message("Data file clinical_cases_byISO_medianCI.csv does not exist in the intermediate_table folder.")
}else{
    clinical_cases_byISO_medianCI <- read.csv(paste0(output_final_path, "/intermediate_table/clinical_cases_byISO_medianCI.csv"))
}

n_tested_byISO <- clinical_cases_byISO_medianCI %>%
    mutate(n_rdt_tested_lb = case_when(confirmation_lens == "global-estimate" | confirmation_lens == "no-estimate" ~ 0,
                                       confirmation_lens == "district-estimate" ~ cumu_clinical_cases_lb * prop_sus_rdt_tested)) %>%
    mutate(n_rdt_tested_median = case_when(confirmation_lens == "global-estimate" | confirmation_lens == "no-estimate" ~ 0,
                                           confirmation_lens == "district-estimate" ~ cumu_clinical_cases_median * prop_sus_rdt_tested)) %>%
    mutate(n_rdt_tested_ub = case_when(confirmation_lens == "global-estimate" | confirmation_lens == "no-estimate" ~ 0,
                                       confirmation_lens == "district-estimate" ~ cumu_clinical_cases_ub * prop_sus_rdt_tested)) %>%
    mutate(n_culture_tested_lb = case_when(confirmation_lens == "no-estimate" ~ 0,
                                           confirmation_lens == "global-estimate" ~ cumu_clinical_cases_lb * prop_sus_cultured,
                                           confirmation_lens == "district-estimate" ~ cumu_clinical_cases_lb * rdt_pos_rate * prop_rdt_positive_cultured)) %>%
    mutate(n_culture_tested_median = case_when(confirmation_lens == "no-estimate" ~ 0,
                                               confirmation_lens == "global-estimate" ~ cumu_clinical_cases_median * prop_sus_cultured,
                                               confirmation_lens == "district-estimate" ~ cumu_clinical_cases_median * rdt_pos_rate * prop_rdt_positive_cultured)) %>%
    mutate(n_culture_tested_ub = case_when(confirmation_lens == "no-estimate" ~ 0,
                                           confirmation_lens == "global-estimate" ~ cumu_clinical_cases_ub * prop_sus_cultured,
                                           confirmation_lens == "district-estimate" ~ cumu_clinical_cases_ub * rdt_pos_rate * prop_rdt_positive_cultured))
write.csv(n_tested_byISO, paste0(cost_path, "n_tested_byISO.csv"), row.names = F)


    
    
