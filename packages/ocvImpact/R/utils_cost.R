## Functions used for calculate cost metrics.

# Function to get clinical cases with median and CI (across years) from one country's target_table in intermediate folder
get_clinical_cases_medianCI <- function(df_tt){

    df_clinical <- df_tt %>%
        mutate(clinical_cases = confirmed_incidence_rate / confirmation_rate * pop_model) %>%
        group_by(ISO, confirmation_lens, admin_level, threshold, run_id, year) %>%
        summarize(clinical_cases = sum(clinical_cases)) %>%
        group_by(ISO, confirmation_lens, admin_level, threshold, run_id) %>% 
        mutate(cumu_clinical_cases = cumsum(clinical_cases)) %>%
        filter(year == 2035) %>% dplyr::select(-year) %>%
        # get median and CI of cumu_clinical_cases
        group_by(ISO, confirmation_lens, admin_level, threshold) %>%
        summarize(cumu_clinical_cases_lb = quantile(cumu_clinical_cases, 0.025, na.rm = T),
                  cumu_clinical_cases_median = quantile(cumu_clinical_cases, 0.5, na.rm = T),
                  cumu_clinical_cases_ub = quantile(cumu_clinical_cases, 0.975, na.rm = T))
  
    return(df_clinical)
}

# Function to get total clinical cases by run (without summarizing median and CI)
get_clinical_cases <- function(df_tt){

    df_clinical <- df_tt %>%
        mutate(clinical_cases = confirmed_incidence_rate / confirmation_rate * pop_model) %>%
        group_by(ISO, confirmation_lens, admin_level, threshold, run_id, year) %>%
        summarize(clinical_cases = sum(clinical_cases)) %>%
        group_by(ISO, confirmation_lens, admin_level, threshold, run_id) %>% 
        mutate(cumu_clinical_cases = cumsum(clinical_cases)) %>%
        filter(year == 2035) %>% dplyr::select(-year, -clinical_cases)
  
    return(df_clinical)
}


# Function to get cost and cost effectiveness from target_table of one country
calc_costeff_ocv_test <- function(df_tt){

    df_clinical <- get_clinical_cases(df_tt)
    
    df_fvp <- df_tt %>%
        group_by(confirmation_lens, admin_level, threshold, year, run_id) %>%
        summarize(fvp = sum(actual_fvp), ac = sum(true_averted_cases)) %>%
        group_by(confirmation_lens, admin_level, threshold, run_id) %>%
        mutate(cumu_fvp = cumsum(fvp), cumu_ac = cumsum(ac)) %>%
        filter(year == 2035) %>%
        dplyr::select(-year, -fvp, -ac)

    df_cost <- df_clinical %>%
        left_join(df_fvp, by = c("confirmation_lens", "admin_level", "threshold", "run_id")) %>%
        mutate(cumu_ac = ifelse(cumu_ac < 1, 1, cumu_ac)) %>% # when the ac < 0, replace it with "one averted case"
        mutate(ocv_cost = cumu_fvp * dose_regimen * cost_per_ocv) %>%
        mutate(rdt_cost = case_when(confirmation_lens == "global-estimate" | confirmation_lens == "no-estimate" ~ 0,
                                    confirmation_lens == "district-estimate" ~ cumu_clinical_cases * prop_sus_rdt_tested * cost_per_rdt)) %>%
        mutate(culture_cost = case_when(confirmation_lens == "no-estimate" ~ 0,
                                        confirmation_lens == "global-estimate" ~ cumu_clinical_cases * prop_sus_cultured * cost_per_culture,
                                        confirmation_lens == "district-estimate" ~ cumu_clinical_cases * rdt_pos_rate * prop_rdt_positive_cultured * cost_per_culture)) %>%
        mutate(test_cost = rdt_cost + culture_cost) %>%
        mutate(costeff_ocv_only = ocv_cost / cumu_ac) %>%
        mutate(costeff_ocv_test = (ocv_cost + test_cost) / cumu_ac) %>% 
        group_by(ISO, confirmation_lens, admin_level, threshold) %>%
        summarize(costeff_ocv_test_lb = quantile(costeff_ocv_test, 0.025, na.rm = T),
                  costeff_ocv_test_median = quantile(costeff_ocv_test, 0.5, na.rm = T),
                  costeff_ocv_test_ub = quantile(costeff_ocv_test, 0.975, na.rm = T)) 
        
        return(df_cost)

}
