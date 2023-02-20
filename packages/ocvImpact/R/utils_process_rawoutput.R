#' @name make_target_table
#' @title make_target_table
#' @description combine all raw outputs from selected countries to make a combined target table
#' @param countries vector of ISO country code(s)
#' @param thresholds selected thresholds, e.g. all thresholds: thresholds = c(0.001, 2e-04, 1e-04)
#' @param surveillance_scenarios selected scenarios, e.g. all scenarios: surveillance_scenarios = c("no-estimate", "global-estimate", "district-estimate")
#' @param output_raw_path path to the folders where raw output dataframes are stored
#' @return a combined target_table for selected countries (with true cases number of different scenarios)
make_target_table <- function(countries, 
                              thresholds,
                              surveillance_scenarios = c("no-estimate", "global-estimate", "district-estimate"),
                              output_raw_path = "output_raw/202302_survms"){
  
    # loop through all the countries 
    for(country in countries){
      
      # make an empty table to store all output from one single country
      df_country <- as.data.frame(matrix(NA, ncol = 21))
      colnames(df_country) <- c("ISO", "NAME_0", "NAME_1", "NAME_2", "confirmation_lens", "admin_level", "threshold", "year","run_id", "true_confirm_rate", "true_incidence_rate",
                                "confirmation_rate", "confirmed_incidence_rate", "pop_model", "pop_prop", "latest_target_year", "is_target",
                                "actual_prop_vaccinated", "actual_fvp", "campaign_default_true_case", "no_vaccination_true_case")

      ## first, get the names of the files (in format: incid_FALSE_outbk_FALSE_0.001_district-estimate_BEN_rc_admin1)
      prefix <- "incid_FALSE_outbk_FALSE_"
      
      for(admin_suffix in c("_rc_admin1", "_rc_admin2")){
        for(threshold in thresholds){
          for(surveillance_scenario in surveillance_scenarios){
            fn_cd <- paste0(prefix, threshold, "_", surveillance_scenario, "_", country, admin_suffix, ".csv")
            fn_novax <- paste0(prefix, "0.001_district-estimate_", country, admin_suffix, ".csv")
            dir_cd <- paste0(output_raw_path, "/campaign-default/")
            dir_novax <- paste0(output_raw_path, "/no-vaccination/")

            # read in raw output 
            df_cd <- read.csv(paste0(dir_cd, fn_cd))
            df_novax <- read.csv(paste0(dir_novax, fn_novax))
       
            # add useful columns and rename columns
            if(!"NAME_2" %in% names(df_cd)){df_cd <- df_cd %>% mutate(NAME_2 = NA)}
            if(!"NAME_2" %in% names(df_novax)){df_novax <- df_novax %>% mutate(NAME_2 = NA)}
            df_cd <- df_cd %>% 
              mutate(threshold = threshold, admin_level = str_remove(admin_suffix, "_rc_")) %>%
              rename(campaign_default_true_case = true_case)
            df_novax <- df_novax %>% 
              mutate(threshold = threshold, admin_level = str_remove(admin_suffix, "_rc_")) %>%
              rename(no_vaccination_true_case = true_case) %>%
              dplyr::select(ISO, admin_level, NAME_1, NAME_2, threshold, year, run_id, no_vaccination_true_case)

            # join df_cd and df_novax together
            df_cd_novax <- df_cd %>% 
              left_join(df_novax,  by = c("ISO", "admin_level", "NAME_1", "NAME_2", "threshold", "year", "run_id"))
            
            df_country <- rbind(df_country, df_cd_novax) %>% filter(ISO != "NA")
          }
        }
      }

      if(which(countries == country) == 1){df_alloutput <- df_country}
      if(which(countries == country) != 1){df_alloutput <- rbind(df_alloutput, df_country)}

    }

    return(df_alloutput)

}



#' @export
#' @name get_tp_eff_median
#' @title get_tp_eff_median
#' @description summarize median targeted pop and median efficiency of each country
#' @param rc target table (can be target table of one or multiple countries)
get_dose_eff_median <- function(rc){
  
  df_tp_cumu <- rc %>% 
    filter(is_target == 1) %>%
    group_by(ISO, year, run_id, threshold, confirmation_lens, admin_level) %>%
    summarize(target_pop = sum(actual_fvp)) %>%
    group_by(ISO, run_id, threshold, confirmation_lens, admin_level) %>%
    mutate(target_pop_cumu = cumsum(target_pop)) %>% # get sum tp of all years
    dplyr::select(-target_pop) %>%
    group_by(ISO, year, threshold, confirmation_lens, admin_level) %>%
    summarize(tp_cumu_median = median(target_pop_cumu, na.rm = TRUE)) %>% # only keep median for each ISO
    dplyr::select(year, threshold, confirmation_lens, tp_cumu_median, admin_level) %>%
    group_by(ISO, threshold, confirmation_lens, admin_level) %>%
    summarize(tp_cumu_median = max(tp_cumu_median))
  
  # averted true cases per 1000 fvp 
  df_eff <- rc %>% 
    group_by(ISO, year, run_id, threshold, confirmation_lens, admin_level) %>%
    summarize(true_ac = sum(true_averted_cases)) %>%
    group_by(ISO, run_id, threshold, confirmation_lens, admin_level) %>%
    mutate(true_ac_cumu = cumsum(true_ac)) %>%
    dplyr::select(-true_ac) %>%
    # filter(year == max_year)
    group_by(ISO, run_id, threshold, confirmation_lens, admin_level) %>%
    summarize(true_ac_cumu = max(true_ac_cumu)) %>%
    right_join(
      rc %>% 
        group_by(ISO, year, run_id, threshold, confirmation_lens, admin_level) %>%
        summarize(fvp = sum(actual_fvp)) %>%
        group_by(ISO, run_id, threshold, confirmation_lens, admin_level) %>%
        mutate(fvp_cumu = cumsum(fvp)) %>%
        dplyr::select(-fvp) %>%
        group_by(ISO, run_id, threshold, confirmation_lens, admin_level) %>%
        summarize(fvp_cumu = max(fvp_cumu)) ,
      by = c("ISO", "run_id", "threshold", "confirmation_lens", "admin_level")
    )%>%
    # filter out all those fvp == 0 rows to get rid of Inf
    filter(fvp_cumu != 0) %>%
    # calculate efficiency
    mutate(efficiency = true_ac_cumu / fvp_cumu * 1000) %>%
    group_by(ISO, threshold, confirmation_lens, admin_level) %>%
    summarize(efficiency_median = median(efficiency, na.rm = TRUE)) %>% # keep only median for each country
    dplyr::select(ISO, threshold, confirmation_lens,admin_level, efficiency_median)
  
  
  # join two tables 
  df_result <- df_tp_cumu %>%
    left_join(df_eff, by = c("ISO", "threshold", "confirmation_lens", "admin_level")) %>%
    # convert unit of tp to "million"
    mutate(doses = tp_cumu_median/1000000 * 2) %>%
    rename(efficiency = efficiency_median) 
  
  return(df_result)
}

#' @name make_eff_table
#' @title make_eff_table
#' @description make a table summarizing (total) target population and ocv efficiency (at country level)
#' @param rc target table of one country at one admin level
#' @return a table summarizing target population and ocv efficiency for one country at one admin_level (either admin1 or admin2)
#' @export
#' @include
make_eff_table_median <- function(rc){
  
  country <- rc$ISO[1]

    df_tp_cumu <- rc %>% 
      #filter(is_target == 1) %>% (need to include 0 dose even for median calculation)
      group_by(year, run_id, threshold, confirmation_lens) %>%
      #summarize(target_pop = sum(pop_model)) %>% # change this to "targeted pop", which is fvp
      summarize(target_pop = sum(actual_fvp)) %>%
      group_by(run_id, threshold, confirmation_lens) %>%
      mutate(target_pop_cumu = cumsum(target_pop)) %>%
      dplyr::select(-target_pop) %>%
      group_by(year, threshold, confirmation_lens) %>%
      summarize(tp_cumu_lb = quantile(target_pop_cumu, 0.025, na.rm = TRUE), # lower 95%CI
                tp_cumu_median = median(target_pop_cumu, na.rm = TRUE),
                tp_cumu_ub = quantile(target_pop_cumu, 0.975, na.rm = TRUE)) %>% # upper 95%CI
      mutate(ISO = country) %>%
      dplyr::select(ISO, year, threshold, confirmation_lens, tp_cumu_lb, tp_cumu_median, tp_cumu_ub) %>%
      group_by(ISO, threshold, confirmation_lens) %>%
      summarize(tp_cumu_lb = max(tp_cumu_lb, na.rm = T),
                tp_cumu_median = max(tp_cumu_median, na.rm = T),
                tp_cumu_ub = max(tp_cumu_ub, na.rm = T))
  
    # efficiency: averted true cases per 1000 fvp 
    df_eff <- rc %>% 
      group_by(year, run_id, threshold, confirmation_lens) %>%
      summarize(true_ac = sum(true_averted_cases)) %>%
      group_by(run_id, threshold, confirmation_lens) %>%
      mutate(true_ac_cumu = cumsum(true_ac)) %>%
      dplyr::select(-true_ac) %>%
      # filter(year == max_year)
      group_by(run_id, threshold, confirmation_lens) %>%
      summarize(true_ac_cumu = max(true_ac_cumu)) %>%
      right_join(
        rc %>% 
          group_by(year, run_id, threshold, confirmation_lens) %>%
          summarize(fvp = sum(actual_fvp)) %>%
          group_by(run_id, threshold, confirmation_lens) %>%
          mutate(fvp_cumu = cumsum(fvp)) %>%
          dplyr::select(-fvp) %>%
          group_by(run_id, threshold, confirmation_lens) %>%
          summarize(fvp_cumu = max(fvp_cumu)) ,
        by = c("run_id", "threshold", "confirmation_lens")
      )%>%
      # added 10/11:
      # filter out all those fvp == 0 rows to get rid of Inf
      filter(fvp_cumu != 0) %>%
      # calculate efficiency
      mutate(efficiency = true_ac_cumu / fvp_cumu * 1000) %>%
      group_by(threshold, confirmation_lens) %>%
      summarize(efficiency_lb = quantile(efficiency, 0.025, na.rm = TRUE), # lower 95%CI
                efficiency_median = median(efficiency, na.rm = TRUE),
                efficiency_ub = quantile(efficiency, 0.975, na.rm = TRUE)) %>% # upper 95%CI
      mutate(ISO = country) %>%
      dplyr::select(ISO, threshold, confirmation_lens, efficiency_lb, efficiency_median, efficiency_ub)
    
    # averted cases:
    df_ac <- rc %>% 
      group_by(year, run_id, threshold, confirmation_lens) %>%
      summarize(true_ac = sum(true_averted_cases)) %>%
      group_by(run_id, threshold, confirmation_lens) %>%
      mutate(true_ac_cumu = cumsum(true_ac)) %>%
      dplyr::select(-true_ac) %>%
      # filter(year == max_year)
      group_by(run_id, threshold, confirmation_lens) %>%
      summarize(true_ac_cumu = max(true_ac_cumu)) %>%
      # calculate median
      group_by(threshold, confirmation_lens) %>%      
      summarize(ac_lb = quantile(true_ac_cumu, 0.025, na.rm = TRUE), # lower 95%CI
                ac_median = median(true_ac_cumu, na.rm = TRUE),
                ac_ub = quantile(true_ac_cumu, 0.975, na.rm = TRUE)) %>% # upper 95%CI
      mutate(ISO = country) %>%
      dplyr::select(ISO, threshold, confirmation_lens, ac_lb, ac_median, ac_ub)
    

    # join two tables and round up numbers
    df_result <- df_tp_cumu %>%
      full_join(df_eff, by = c("ISO", "threshold", "confirmation_lens")) %>%
      full_join(df_ac, by = c("ISO", "threshold", "confirmation_lens"))
      # convert unit of tp to "million"
      #mutate(tp_cumu_lb = round(tp_cumu_lb/1000000, 2),
      #       tp_cumu_median = round(tp_cumu_median/1000000, 2),
      #       tp_cumu_ub = round(tp_cumu_ub/1000000, 2),
      #       efficiency_lb = round(efficiency_lb, 2),
      #       efficiency_median = round(efficiency_median, 2),
      #       efficiency_ub = round(efficiency_ub, 2))
    

    # format df_result to "median(lb-ub)" format
    #df_result <- df_result %>%
    #  mutate(tp_cumu = paste0(tp_cumu_median, " (", tp_cumu_lb, "-", tp_cumu_ub, ")"),
    #         efficiency = paste0(efficiency_median, " (", efficiency_lb, "-", efficiency_ub, ")"),
    #         ac = paste0(ac_median, " (", efficiency_lb, "-", efficiency_ub, ")")) %>%
    #  select(ISO, threshold, confirmation_lens, tp_cumu, efficiency, ac)
  
  return(df_result)
}


#' @name combine_eff_table
#' @title combine_eff_table
#' @description combine tables from multiple countries and both admin levels
#' @param countries a vector contains all the countries included in the diagnostic report 
#' @param df_target target_table generated and saved in cache
#' @param admin_level admin level(s) selected for the diagnostic report
#' @return a table summarizing target population and ocv efficiency for one or more countries (at country level)
#' @export
#' @include
combine_eff_table_median <- function(countries, 
                              df_target,
                              admin_level){
    
    for(country in countries){
      
      if(admin_level %in% c("both", "admin1")){
        df_admin1 <- df_target %>% filter(admin_level == "admin1" & ISO == country)
        eff_table_admin1 <- make_eff_table_median(df_admin1) %>% mutate(admin_level = "admin1")
      }
      if(admin_level %in% c("both", "admin2")){
        df_admin2 <- df_target %>% filter(admin_level == "admin2" & ISO == country)
        eff_table_admin2 <- make_eff_table_median(df_admin2) %>% mutate(admin_level = "admin2")
      }

      # combine results
      if(admin_level == "admin1"){eff_table <- eff_table_admin1}
      if(admin_level == "admin2"){eff_table <- eff_table_admin2}
      if(admin_level == "both"){eff_table <- rbind(eff_table_admin1, eff_table_admin2)}
      
      # if this is the first country
      if(which(countries == country) == 1){df_eff_comb <- eff_table}
      if(which(countries == country) != 1){df_eff_comb <- rbind(df_eff_comb, eff_table)}

    }

    return(df_eff_comb)
}

#' @name make_eff_table_allISOs
#' @title make_eff_table_allISOs
#' @description make a table summarizing (total) target population and ocv efficiency, combining all the countries all at once
#' @param rc target table all countries at one admin level
#' @return a table summarizing target population and ocv efficiency for all countries at one admin_level (either admin1 or admin2)
#' @export
#' @include
make_eff_table_allISOs <- function(rc){
  
    df_tp_cumu <- rc %>% 
      filter(is_target == 1) %>%
      group_by(year, run_id, threshold, confirmation_lens) %>%
      summarize(target_pop = sum(actual_fvp)) %>%
      group_by(run_id, threshold, confirmation_lens) %>%
      mutate(target_pop_cumu = cumsum(target_pop)) %>%
      dplyr::select(-target_pop) %>%
      group_by(year, threshold, confirmation_lens) %>%
      summarize(tp_cumu_lb = quantile(target_pop_cumu, 0.025, na.rm = TRUE), # lower 95%CI
                tp_cumu_median = median(target_pop_cumu, na.rm = TRUE),
                tp_cumu_ub = quantile(target_pop_cumu, 0.975, na.rm = TRUE)) %>% # upper 95%CI
      dplyr::select(year, threshold, confirmation_lens, tp_cumu_lb, tp_cumu_median, tp_cumu_ub) %>%
      group_by(threshold, confirmation_lens) %>%
      summarize(tp_cumu_lb = max(tp_cumu_lb),
                tp_cumu_median = max(tp_cumu_median),
                tp_cumu_ub = max(tp_cumu_ub))
  
    # averted true cases per 1000 fvp 
    df_eff <- rc %>% 
      group_by(year, run_id, threshold, confirmation_lens) %>%
      summarize(true_ac = sum(true_averted_cases)) %>%
      group_by(run_id, threshold, confirmation_lens) %>%
      mutate(true_ac_cumu = cumsum(true_ac)) %>%
      dplyr::select(-true_ac) %>%
      # filter(year == max_year)
      group_by(run_id, threshold, confirmation_lens) %>%
      summarize(true_ac_cumu = max(true_ac_cumu)) %>%
      right_join(
        rc %>% 
          group_by(year, run_id, threshold, confirmation_lens) %>%
          summarize(fvp = sum(actual_fvp)) %>%
          group_by(run_id, threshold, confirmation_lens) %>%
          mutate(fvp_cumu = cumsum(fvp)) %>%
          dplyr::select(-fvp) %>%
          group_by(run_id, threshold, confirmation_lens) %>%
          summarize(fvp_cumu = max(fvp_cumu)) ,
        by = c("run_id", "threshold", "confirmation_lens")
      )%>%
      # added on 10/11:
      # filter out all those fvp == 0 rows to get rid of Inf
      filter(fvp_cumu != 0) %>%
      # calculate efficiency
      mutate(efficiency = true_ac_cumu / fvp_cumu * 1000) %>%
      group_by(threshold, confirmation_lens) %>%
      summarize(efficiency_lb = quantile(efficiency, 0.025, na.rm = TRUE), # lower 95%CI
                efficiency_median = median(efficiency, na.rm = TRUE),
                efficiency_ub = quantile(efficiency, 0.975, na.rm = TRUE)) %>% # upper 95%CI
      dplyr::select(threshold, confirmation_lens, efficiency_lb, efficiency_median, efficiency_ub)
    
    
    # join two tables and round up numbers
    df_result <- df_tp_cumu %>%
      left_join(df_eff, by = c("threshold", "confirmation_lens")) %>%
      # convert unit of tp to "million"
      mutate(tp_cumu_lb = round(tp_cumu_lb/1000000, 2),
             tp_cumu_median = round(tp_cumu_median/1000000, 2),
             tp_cumu_ub = round(tp_cumu_ub/1000000, 2),
             efficiency_lb = round(efficiency_lb, 2),
             efficiency_median = round(efficiency_median, 2),
             efficiency_ub = round(efficiency_ub, 2))
    
    # format df_result to "median(lb-ub)" format
    df_result <- df_result %>%
      mutate(tp_cumu = paste0(tp_cumu_median, " (", tp_cumu_lb, "-", tp_cumu_ub, ")"),
             efficiency = paste0(efficiency_median, " (", efficiency_lb, "-", efficiency_ub, ")")) %>%
      select(threshold, confirmation_lens, tp_cumu, efficiency)
  
  return(df_result)
}

#' @name combine_eff_table_allISOs
#' @title combine_eff_table_allISOs
#' @description combine tables from different admin level(s)
#' @param countries a vector contains all the countries included in the diagnostic report 
#' @param df_target target_table generated and saved in cache
#' @param admin_level admin level(s) selected for the diagnostic report
#' @return a table summarizing target population and ocv efficiency for one or more countries (at country level)
#' @export
#' @include
combine_eff_table_allISOs <- function(df_target,
                                      admin_level){
    
      if(admin_level %in% c("both", "admin1")){
        df_admin1 <- df_target %>% filter(admin_level == "admin1")
        eff_table_admin1 <- make_eff_table_allISOs(df_admin1) %>% mutate(admin_level = "admin1")
      }
      if(admin_level %in% c("both", "admin2")){
        df_admin2 <- df_target %>% filter(admin_level == "admin2")
        eff_table_admin2 <- make_eff_table_allISOs(df_admin2) %>% mutate(admin_level = "admin2")
      }

      # combine results
      if(admin_level == "admin1"){eff_table <- eff_table_admin1}
      if(admin_level == "admin2"){eff_table <- eff_table_admin2}
      if(admin_level == "both"){eff_table <- rbind(eff_table_admin1, eff_table_admin2)}
      
    return(eff_table)
}


## count how many of runs are having no targeting going on for each country
count_targeted_runs <- function(rc,
                                 country){
    df_admin1 <- rc %>% filter(admin_level == "admin1")
    df_admin2 <- rc %>% filter(admin_level == "admin2")

    n_targeted_admin1 <- df_admin1 %>% 
          group_by(year, run_id, threshold, confirmation_lens) %>%
          summarize(fvp = sum(actual_fvp)) %>%
          group_by(run_id, threshold, confirmation_lens) %>%
          mutate(fvp_cumu = cumsum(fvp)) %>%
          dplyr::select(-fvp) %>%
          group_by(run_id, threshold, confirmation_lens) %>%
          summarize(fvp_cumu = max(fvp_cumu)) %>%
          group_by(threshold, confirmation_lens) %>%
          filter(fvp_cumu != 0) %>%
          count() %>%
          mutate(admin_level = "admin1")

    n_targeted_admin2 <- df_admin2 %>% 
          group_by(year, run_id, threshold, confirmation_lens) %>%
          summarize(fvp = sum(actual_fvp)) %>%
          group_by(run_id, threshold, confirmation_lens) %>%
          mutate(fvp_cumu = cumsum(fvp)) %>%
          dplyr::select(-fvp) %>%
          group_by(run_id, threshold, confirmation_lens) %>%
          summarize(fvp_cumu = max(fvp_cumu)) %>%
          group_by(threshold, confirmation_lens) %>%
          filter(fvp_cumu != 0) %>%
          count() %>%
          mutate(admin_level = "admin2")
    
    n_targeted <- rbind(n_targeted_admin1, n_targeted_admin2)
    
    df_n <- expand.grid(ISO = country, threshold = c(0.001, 2e-04, 1e-04), confirmation_lens = c("district-estimate", "global-estimate", "no-estimate"),
                admin_level = c("admin1", "admin2"))
    df_n <- df_n %>% left_join(n_targeted, by = c("threshold", "confirmation_lens", "admin_level")) %>%
                     mutate(n_runs_targeted = ifelse(is.na(n), 0, n)) %>% dplyr::select(-n)
    
    return(df_n)


}


