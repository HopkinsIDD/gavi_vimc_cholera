###################
# Functions
###################

#' @name check_model_setting
#' @title check_model_setting
#' @description Table of model settings & assumptions (ie, which scenario does the diagnostic report represent?)
#' @param rawoutpath (might not need this)
#' @param configpath path of configuration files
#' @param countries vector of all the countries included the diagnostic report
#' @param scenario takes value of either "campaign-default" or "no-vaccination"
#' @param threshold threshold chosen for the diagnostic report, takes one of 1/1000, 1/5000, or 1/10,000
#' @param surveillance_scenario scale of surveillance, takes "no-estimate", "global-estimate", or "district-estimate"
#' @return table with model settings & scenario parameters
#' @export
#' @include
check_model_setting <- function(configpath = "configs/202110gavi-3", 
                                countries,
                                scenario,   
                                surveillance_scenario, 
                                thresholds, 
                                ...){   
  
  # make a table for the output
  df_setting <- matrix(NA, nrow = length(countries), ncol = 13) %>% as.data.frame()
  names(df_setting) <- c("country", "scenario", "surveillance_scenario", "threshold", "num_samples", "vac_coverage", "num_skip_years", "sim_start_year", "vac_start_year", "vac_end_year", "sim_end_year", "incidence_trend", "outbreak_multiplier")
  df_setting$country <- countries

  # check if the desiring config files are available in config folder
  for (threshold in thresholds){

    for(country in countries){
      config_file_path <- paste0(configpath, "/", scenario, "/", surveillance_scenario, "/", threshold)
      config_file_name <- paste0(country, "_", scenario, "_", surveillance_scenario)
      
      if(length(list.files(path = config_file_path, pattern = paste0("^", config_file_name))) == 0){
        message(paste("The configuration file for", country, "under scenario of", scenario, "using", surveillance_scenario, "for confirmation rate and a threshold of", threshold, "doesn't exist in", config_file_path, "."))
      }else{
        config <- yaml::read_yaml(paste0(config_file_path, "/", config_file_name, "_50.yml"))
        df_setting[which(df_setting$country == country),]$scenario <- config$scenario
        df_setting[which(df_setting$country == country),]$num_samples <- as.numeric(config$vacc$num_skip_years) 
        df_setting[which(df_setting$country == country),]$threshold <- as.numeric(config$vacc$vac_incid_threshold)
        df_setting[which(df_setting$country == country),]$vac_coverage <- as.numeric(config$vacc$vac_coverage)
        df_setting[which(df_setting$country == country),]$num_skip_years <- as.numeric(config$vacc$num_skip_years)
        df_setting[which(df_setting$country == country),]$sim_start_year <- as.numeric(config$vacc$sim_start_year)
        df_setting[which(df_setting$country == country),]$vac_start_year <- as.numeric(config$vacc$vac_start_year)
        df_setting[which(df_setting$country == country),]$vac_end_year <- as.numeric(config$vacc$vac_end_year)
        df_setting[which(df_setting$country == country),]$sim_end_year <- as.numeric(config$vacc$sim_end_year)
        df_setting[which(df_setting$country == country),]$surveillance_scenario <- config$surveillance_scenario$surveillance_scenario
        df_setting[which(df_setting$country == country),]$incidence_trend <- config$setting$incidence_rate_trend
        df_setting[which(df_setting$country == country),]$outbreak_multiplier <- config$setting$outbreak_multiplier
      }
    }
    if(!exists("df_setting_final")){
      df_setting_final <- df_setting
    }else{
      df_setting_final <- rbind(df_setting_final, df_setting)
    }

  }
  return(df_setting_final)
}



#' @name check_outputs_availability
#' @title check_outputs_availability
#' @description function to check if all the output files are there 
#' @param rawoutpath path of output_raw in which all the model output files are stored
#' @param configpath path of configuration files
#' @param countries vector of all the countries included the diagnostic report
#' @param scenario takes value of either "campaign-default" or "no-vaccination"
#' @param threshold threshold chosen for the diagnostic report, takes one of 1/1000, 1/5000, or 1/10,000
#' @param surveillance_scenario scale of surveillance, takes "no-estimate", "global-estimate", or "district-estimate"
#' @param admin_level the admin level to check for output
#' @return TRUE or FALSE indicating whether the output files for the specified countries and scenarios exist in output folder, and also return warning messages to check output of which country is missing
#' @export
#' @include
check_outputs_availability <- function(rawoutpath = "output_raw/202110gavi-3",
                                        configpath = "configs/202110gavi-3",
                                        countries,
                                        scenario,
                                        surveillance_scenario,
                                        threshold,
                                        admin_level, 
                                        ...){
  
  df_setting <- check_model_setting(configpath = configpath, countries = countries,
                                scenario = scenario, surveillance_scenario = surveillance_scenario, 
                                threshold = threshold)

  log <- rep(NA, length(countries))
  output_path <- paste0(rawoutpath, "/", scenario)

  for(country in countries){

    incid <- df_setting[which(df_setting$country == country), ]$incidence_trend 
    outbk <- df_setting[which(df_setting$country == country), ]$outbreak_multiplier 
    prefix <- paste0("incid_", incid, "_outbk_", outbk, "_")
    
    filename_admin1 <- paste0(output_path, "/", prefix, threshold, "_", surveillance_scenario, "_", country, "_rc_admin1.csv")
    filename_admin2 <- paste0(output_path, "/", prefix, threshold, "_", surveillance_scenario, "_", country, "_rc_admin2.csv")
    
    if(   ((admin_level == "both")&(!file.exists(filename_admin1) | !file.exists(filename_admin2)))   | 
          (admin_level == "admin1"&!file.exists(filename_admin1))   | 
          (admin_level == "admin2"&!file.exists(filename_admin2))   ){
      message(paste("Either admin1 or admin2 level output is missing for", country, "under scenario of", scenario, "with", surveillance_scenario, "of confirmation rate and threshold of", threshold, "."))
      log[which(countries == country)] <- FALSE
    }else{
      log[which(countries == country)] <- TRUE
    }
  }

  if(FALSE %in% log){
      return(FALSE)
    }else{
      return(TRUE)
    }

}



#' @export
#' @name get_filenames
#' @title get_filenames
#' @description plot the polygon with modeled cases
#' @param cache the cached environment
#' @param surveillance_project_directory surveillance project directory
#' @param pre_country countries specified in the parameters 
#' @param pre_vac_incid_thresholds incidence rate threshold specified in the parameters
#' @param pre_admin_levels admin levels specified in the parameters
#' @param no_vaccination_surveillance_scenario the surveillance scenario specified for the no vaccination simulation 
#' @return cached file names 
get_filenames <- function(cache, surveillance_project_directory, pre_country = NULL, pre_vac_incid_thresholds = NULL, pre_admin_levels = NULL, 
                          no_vaccination_surveillance_scenario){
  
  # Get the set_all_parameters.R to have the context first
  source(file.path(surveillance_project_directory, "scripts/set_all_parameters.R"))

  
  # Get the correct version of country list and threshold list
  if(is.null(pre_country)){
    country_list <- ifelse(pre_country %in% countries, pre_country, countries)  
    rm(pre_country)
  }else{country_list <- pre_country}

  if(is.null(pre_vac_incid_thresholds)){
    threshold_list <- ifelse(as.numeric(pre_vac_incid_thresholds) %in% vac_incid_thresholds, pre_vac_incid_thresholds, vac_incid_thresholds)
    rm(pre_vac_incid_thresholds)
  }else{threshold_list <- pre_vac_incid_thresholds}

  if(!is.null(pre_admin_levels)){vac_admin_level <- pre_admin_levels}
  if(identical(vac_admin_level, "both")) {admin_list <- c("admin1", "admin2")} else {admin_list <- vac_admin_level}
  year_list <- sim_start_year:sim_end_year

  
  # Generate and cache the file names 
  for(i in 1:length(country_list)){
    ## expected cases raster
    ec_filenames_camp <- apply(expand.grid(paste0("output_raw/", runname, "/campaign-default/", "incid"), 
                              incidence_rate_trend, "outbk", outbreak_multiplier, threshold_list, surveillance_scenarios, 
                              country_list[i], "ec", admin_list, paste0(year_list, ".tif")), 1, paste, collapse = "_")
    ec_filenames_novac <- apply(expand.grid(paste0("output_raw/", runname, "/no-vaccination/", "incid"), 
                              incidence_rate_trend, "outbk", outbreak_multiplier, threshold_list, no_vaccination_surveillance_scenario, 
                              country_list[i], "ec", admin_list, paste0(year_list, ".tif")), 1, paste, collapse = "_")
    ecraster_fns <- c(ec_filenames_camp, ec_filenames_novac)
    ecraster_fns <- paste0(surveillance_project_directory, "/", ecraster_fns)
    # if(sum(!file.exists(ecraster_fns)) > 0){
    #   message(paste0("Cannot find the following files: ", ecraster_fns[!file.exists(ecraster_fns)]))
    #   stop(paste("Incomplete model outputs for expected cases rasters for country", country_list[i]))
    # }

    ## expected cases table
    ectable_fns_camp <- apply(expand.grid(paste0("output_raw/", runname, "/campaign-default/", country_list[i]), 
                              "ec", paste0(admin_list, ".csv")), 1, paste, collapse = "_")
    ectable_fns_novac <- apply(expand.grid(paste0("output_raw/", runname, "/no-vaccination/", country_list[i]), 
                              "ec", paste0(admin_list, ".csv")), 1, paste, collapse = "_")                          
    ectable_fns <- c(ectable_fns_camp, ectable_fns_novac)
    ectable_fns <- paste0(surveillance_project_directory, "/", ectable_fns)
    # if(sum(!file.exists(ectable_fns)) > 0){stop(paste("Incomplete model outputs for expected cases table for country", country_list[i]))}

    ## target table
    targettable_fns_camp <- apply(expand.grid(paste0("output_raw/", runname, "/campaign-default/", "incid"), 
                              incidence_rate_trend, "outbk", outbreak_multiplier, threshold_list, surveillance_scenarios, 
                              country_list[i], "rc", paste0(admin_list, ".csv")), 1, paste, collapse = "_")
    targettable_fns_novac <- apply(expand.grid(paste0("output_raw/", runname, "/no-vaccination/", "incid"), 
                              incidence_rate_trend, "outbk", outbreak_multiplier, threshold_list, no_vaccination_surveillance_scenario, 
                              country_list[i], "rc", paste0(admin_list, ".csv")), 1, paste, collapse = "_")
    targettable_fns <- c(targettable_fns_camp, targettable_fns_novac)
    targettable_fns <- paste0(surveillance_project_directory, "/", targettable_fns)
    if(sum(!file.exists(targettable_fns)) > 0){stop(paste("Incomplete model outputs for target tables for country", country_list[i]))}

    ## time table
    timetable_fns_camp <- apply(expand.grid(paste0("output_raw/", runname, "/campaign-default/", "incid"), 
                              incidence_rate_trend, "outbk", outbreak_multiplier, threshold_list, surveillance_scenarios, 
                              country_list[i], "time_table.csv"), 1, paste, collapse = "_")
    timetable_fns_novac <- apply(expand.grid(paste0("output_raw/", runname, "/no-vaccination/", "incid"), 
                              incidence_rate_trend, "outbk", outbreak_multiplier, threshold_list, no_vaccination_surveillance_scenario, 
                              country_list[i], "time_table.csv"), 1, paste, collapse = "_")
    timetable_fns <- c(timetable_fns_camp, timetable_fns_novac)
    timetable_fns <- paste0(surveillance_project_directory, "/", timetable_fns)
    if(sum(!file.exists(timetable_fns)) > 0){stop(paste("Incomplete model outputs for time tables for country", country_list[i]))}

    ## true confirmation rate raster
    truecomfirmrate_fns <- paste0("intermediate_raster/", country_list[i], "_trueconfirmrate_admin2.tif")  
    truecomfirmrate_fns <- paste0(surveillance_project_directory, "/", truecomfirmrate_fns)
    # if(sum(!file.exists(truecomfirmrate_fns)) > 0){stop(paste("Incomplete model outputs for true confirmation rate raster for country", country_list[i]))}                        

    ## cache all the filenames
    cache[[country_list[i]]]$cases_raster_fns <- ecraster_fns
    cache[[country_list[i]]]$cases_table_fns <- ectable_fns
    cache[[country_list[i]]]$target_table_fns <- targettable_fns
    cache[[country_list[i]]]$time_table_fns <- timetable_fns
    cache[[country_list[i]]]$confirm_rate_raster_fns <- truecomfirmrate_fns
  }
}



#' @export
#' @name combine_output
#' @title combine_output
#' @description combine the output from different countries 
#' @param cache the cached environment
#' @param surveillance_project_directory surveillance project directory
#' @param output_to_combine output to combine
#' @param save_final_output whether or not to save the combined results in the final output folder 
#' @return cached combined results
combine_output <- function(cache, surveillance_project_directory, output_to_combine, save_final_output = TRUE){
  
  if("target_table" %in% output_to_combine){
    # initialize
    target_table <- readr::read_csv(cache[[names(cache)[1]]]$target_table_fns[1])
    target_table <- target_table[c(), ] %>% mutate(NAME_2 = as.character(), incid_trend = as.logical(), outbk_trend = as.logical(), 
      threshold = as.numeric(), surveillance_scenario = as.character(), admin_level = as.character(), campaign_default_true_case = as.numeric())
    tmp_novac <- target_table

    # loop through countries and scenarios 
    for(country in names(cache)){
      for(gen_scn in c("campaign-default", "no-vaccination")){ #first get all the campaign data in
        for(fn in cache[[country]]$target_table_fns){
          
          general_scenario <- ifelse(grepl("campaign-default", fn), "campaign-default", "no-vaccination")
          if(general_scenario != gen_scn){next}
          fn_in_list <- unlist(strsplit(gsub(paste0(surveillance_project_directory, "/"), "", fn), split = "_"))
          tmp <- readr::read_csv(fn)

          if(general_scenario == "campaign-default"){
            target_table <- rbind(target_table, tmp %>% { if(!"NAME_2" %in% names(tmp)) mutate(., NAME_2 = NA) else . } %>% 
                                                        mutate(incid_trend = as.logical(fn_in_list[3]), outbk_trend = as.logical(fn_in_list[5]), 
                                                              threshold = as.numeric(fn_in_list[6]), surveillance_scenario = fn_in_list[7], admin_level = gsub(".csv", "", fn_in_list[10]), 
                                                              ) %>% 
                                                        rename(campaign_default_true_case = true_case))
          }else{
            tmp_novac <- rbind(tmp_novac, tmp %>% { if(!"NAME_2" %in% names(tmp)) mutate(., NAME_2 = NA) else . } %>% 
                                                  mutate(incid_trend = as.logical(fn_in_list[3]), outbk_trend = as.logical(fn_in_list[5]), 
                                                        threshold = as.numeric(fn_in_list[6]), surveillance_scenario = fn_in_list[7], admin_level = gsub(".csv", "", fn_in_list[10])
                                                        ) %>% 
                                                  rename(no_vaccination_true_case = true_case) %>% 
                                                  dplyr::select(ISO, admin_level, NAME_1, NAME_2, incid_trend, outbk_trend, threshold, year, run_id, no_vaccination_true_case))
          }
                                                         
        }
      }
    }
    target_table <- dplyr::left_join(target_table, tmp_novac, by = c("ISO", "admin_level", "NAME_1", "NAME_2", "incid_trend", "outbk_trend", "threshold", "year", "run_id"))
    rm(tmp_novac)
    
    # save the table in the final output folder 
    if(save_final_output){
      runname <- unlist(strsplit(gsub(paste0(surveillance_project_directory, "/"), "", fn), "/"))[2]
      readr::write_csv(target_table, paste0("output_final/", runname, "/target_table.csv"))
    }
    cache$target_table <- target_table
    rm(target_table)
  }


  if("time_table" %in% output_to_combine){
    # initialize
    time_table <- readr::read_csv(cache[[names(cache)[1]]]$time_table_fns[1])
    time_table <- time_table[c(), ] %>% mutate( general_scenario = as.character(), incid_trend = as.logical(), outbk_trend = as.logical(), 
                                                threshold = as.numeric(), surveillance_scenario = as.character(), country = as.character())
    # loop through countries and scenarios 
    for(country in names(cache)){
      for(fn in cache[[country]]$time_table_fns){
        general_scenario <- ifelse(grepl("campaign-default", fn), "campaign-default", "no-vaccination")
        fn_in_list <- unlist(strsplit(gsub(paste0(surveillance_project_directory, "/"), "", fn), split = "_"))
        tmp <- readr::read_csv(fn)
        time_table <- rbind(time_table, tmp %>% mutate( general_scenario = general_scenario, incid_trend = as.logical(fn_in_list[3]), outbk_trend = as.logical(fn_in_list[5]), 
                                                        threshold = as.numeric(fn_in_list[6]), surveillance_scenario = fn_in_list[7], country = country))
      }
    }
    # save the table in the final output folder 
    time_table <- time_table %>% mutate_at(vars(load_baseline_incidence, update_targets_list, update_save_vac_raster, 
                                                update_save_sus_raster,create_expectedCases,add_new_row_to_target_list), funs(round(., 2)))

    if(save_final_output){
      runname <- unlist(strsplit(gsub(paste0(surveillance_project_directory, "/"), "", fn), "/"))[2]
      readr::write_csv(time_table, paste0("output_final/", runname, "/time_table.csv"))
    }
    cache$time_table <- time_table
    rm(time_table)
  }

}



#' @export
#' @name plot_cases
#' @title plot_cases
#' @description plot the cases and return the ggplot object
#' @param cache the cached environment
#' @param case_type type of the cases to plot
#' @param threshold threshold to plot  
#' @param cumulative_type either by year (non_cumulative), or cumulated across years, or cumulated within each district
#' @param no_vaccination_surveillance_scenario the surveillance scenario for the no vaccination scenario 
#' @return cached file names 
plot_cases <- function(cache, case_type, threshold, cumulative_type, no_vaccination_surveillance_scenario){
  
  # Aviod the same name
  chosen_threshold <- threshold

  # Restructure of the target table 
  if(grepl("true", case_type) | grepl("clinical", case_type) | grepl("confirmed", case_type)){
    if(!"target_table_long" %in% names(cache)){
      target_table <- cache$target_table %>% 
        rename(campaign_default = campaign_default_true_case) %>%
        rename(no_vaccination = no_vaccination_true_case) 
      target_table <- target_table %>% 
        tidyr::gather("general_scenario", "true_case", match(c("campaign_default", "no_vaccination"), names(target_table))) %>%
        mutate(clinical_case = true_case / true_confirm_rate) %>% 
        mutate(confirmed_case = clinical_case * confirmation_rate)
      target_table <- target_table %>% 
        mutate(is_target = case_when(general_scenario == "no_vaccination" ~ 0, T ~ is_target)) %>% 
        mutate(actual_prop_vaccinated = case_when(general_scenario == "no_vaccination" ~ 0, T ~ actual_prop_vaccinated)) %>% 
        mutate(actual_fvp = case_when(general_scenario == "no_vaccination" ~ 0, T ~ actual_fvp))
      cache$target_table_long <- target_table
    }else{
      target_table <- cache$target_table_long
    }
    
  }else{
    if(!"target_table_wide" %in% names(cache)){
      target_table <- cache$target_table %>%
        mutate(general_scenario = "both", 
              campaign_default_clinical_case = campaign_default_true_case / true_confirm_rate, 
              no_vaccination_clinical_case = no_vaccination_true_case / true_confirm_rate, 
              campaign_default_confirmed_case = campaign_default_clinical_case * confirmation_rate, 
              no_vaccination_confirmed_case = no_vaccination_clinical_case * confirmation_rate) %>% 
        mutate(averted_true_case = no_vaccination_true_case - campaign_default_true_case, 
              averted_clinical_case = no_vaccination_clinical_case - campaign_default_clinical_case, 
              averted_confirmed_case = no_vaccination_confirmed_case - campaign_default_confirmed_case) 
      cache$target_table_wide <- target_table
    }else{
      target_table <- cache$target_table_wide
    }
  }


  # Add the mean across simulations (sum across districts or years) 
  if(all(c("campaign_default_true_case", "no_vaccination_true_case") %in% names(target_table))){
    if(grepl("non", cumulative_type) | grepl("country", cumulative_type)){
      target_table_total_country <- target_table %>% 
        group_by(ISO, admin_level, incid_trend, outbk_trend, general_scenario, surveillance_scenario, threshold, year, run_id) %>% 
        summarise(campaign_default_true_case = sum(campaign_default_true_case), 
                  no_vaccination_true_case = sum(no_vaccination_true_case), 
                  campaign_default_clinical_case = sum(campaign_default_clinical_case), 
                  no_vaccination_clinical_case = sum(no_vaccination_clinical_case), 
                  campaign_default_confirmed_case = sum(campaign_default_confirmed_case), 
                  no_vaccination_confirmed_case = sum(no_vaccination_confirmed_case), 
                  averted_true_case = sum(averted_true_case), 
                  averted_clinical_case = sum(averted_clinical_case), 
                  averted_confirmed_case = sum(averted_confirmed_case), 
                  dose = sum(actual_fvp*2)) %>%
        mutate(run_id = paste0("total_country", run_id))
      target_table_mean_country <- target_table_total_country %>%
        group_by(ISO, admin_level, incid_trend, outbk_trend, general_scenario, surveillance_scenario, threshold, year) %>% 
        summarise(campaign_default_true_case = mean(campaign_default_true_case), 
                  no_vaccination_true_case = mean(no_vaccination_true_case), 
                  campaign_default_clinical_case = mean(campaign_default_clinical_case), 
                  no_vaccination_clinical_case = mean(no_vaccination_clinical_case), 
                  campaign_default_confirmed_case = mean(campaign_default_confirmed_case), 
                  no_vaccination_confirmed_case = mean(no_vaccination_confirmed_case), 
                  averted_true_case = mean(averted_true_case), 
                  averted_clinical_case = mean(averted_clinical_case), 
                  averted_confirmed_case = mean(averted_confirmed_case), 
                  dose = mean(dose)) %>%
        mutate(run_id = "mean_country")
      target_table_country <- plyr::rbind.fill(target_table_total_country, target_table_mean_country)
      # target_table <- plyr::rbind.fill(target_table, target_table_total_country, target_table_mean_country)
      # rm(target_table_total_country, target_table_mean_country)
    }else{ ### might be only for the district level one
      target_table_total_district_across_year <- target_table %>%
        group_by(ISO, admin_level, NAME_1, NAME_2, incid_trend, outbk_trend, general_scenario, surveillance_scenario, threshold, run_id) %>% 
        summarise(campaign_default_true_case = sum(campaign_default_true_case), 
                  no_vaccination_true_case = sum(no_vaccination_true_case), 
                  campaign_default_clinical_case = sum(campaign_default_clinical_case), 
                  no_vaccination_clinical_case = sum(no_vaccination_clinical_case), 
                  campaign_default_confirmed_case = sum(campaign_default_confirmed_case), 
                  no_vaccination_confirmed_case = sum(no_vaccination_confirmed_case), 
                  averted_true_case = sum(averted_true_case), 
                  averted_clinical_case = sum(averted_clinical_case), 
                  averted_confirmed_case = sum(averted_confirmed_case), 
                  dose = sum(actual_fvp*2)) %>%
        mutate(run_id = "total_district")
      target_table_mean_district_across_year <- target_table_total_district_across_year %>%
        group_by(ISO, admin_level, NAME_1, NAME_2, incid_trend, outbk_trend, general_scenario, surveillance_scenario, threshold) %>% 
        summarise(campaign_default_true_case = mean(campaign_default_true_case), 
                  no_vaccination_true_case = mean(no_vaccination_true_case), 
                  campaign_default_clinical_case = mean(campaign_default_clinical_case), 
                  no_vaccination_clinical_case = mean(no_vaccination_clinical_case), 
                  campaign_default_confirmed_case = mean(campaign_default_confirmed_case), 
                  no_vaccination_confirmed_case = mean(no_vaccination_confirmed_case), 
                  averted_true_case = mean(averted_true_case), 
                  averted_clinical_case = mean(averted_clinical_case), 
                  averted_confirmed_case = mean(averted_confirmed_case), 
                  dose = mean(dose)) %>%
        mutate(run_id = "mean_district")
      target_table_district <- plyr::rbind.fill(target_table_total_district_across_year, target_table_mean_district_across_year)
      # target_table <- plyr::rbind.fill(target_table, target_table_total_district_across_year, target_table_mean_district_across_year)
      # rm(target_table_total_district_across_year, target_table_mean_district_across_year)
    }

  }else{
    if(grepl("non", cumulative_type) | grepl("country", cumulative_type)){
      target_table_total_country <- target_table %>% 
        group_by(ISO, admin_level, incid_trend, outbk_trend, general_scenario, surveillance_scenario, threshold, year, run_id) %>% 
        summarise(true_case = sum(true_case), 
                  clinical_case = sum(clinical_case), 
                  confirmed_case = sum(confirmed_case), 
                  dose = sum(actual_fvp*2)) %>%
        mutate(run_id = paste0("total_country", run_id))
      target_table_mean_country <- target_table_total_country %>%
        group_by(ISO, admin_level, incid_trend, outbk_trend, general_scenario, surveillance_scenario, threshold, year) %>% 
        summarise(true_case = mean(true_case), 
                  clinical_case = mean(clinical_case), 
                  confirmed_case = mean(confirmed_case), 
                  dose = mean(dose)) %>%
        mutate(run_id = "mean_country")
      target_table_country <- plyr::rbind.fill(target_table_total_country, target_table_mean_country)
      # target_table <- plyr::rbind.fill(target_table, target_table_total_country, target_table_mean_country)
      # rm(target_table_total_country, target_table_mean_country)
    }else{
      target_table_total_district_across_year <- target_table %>%
        group_by(ISO, admin_level, NAME_1, NAME_2, incid_trend, outbk_trend, general_scenario, surveillance_scenario, threshold, run_id) %>% 
        summarise(true_case = sum(true_case), 
                  clinical_case = sum(clinical_case), 
                  confirmed_case = sum(confirmed_case), 
                  dose = sum(actual_fvp*2)) %>%
        mutate(run_id = "total_district")
      target_table_mean_district_across_year <- target_table_total_district_across_year %>%
        group_by(ISO, admin_level, NAME_1, NAME_2, incid_trend, outbk_trend, general_scenario, surveillance_scenario, threshold) %>% 
        summarise(true_case = mean(true_case), 
                  clinical_case = mean(clinical_case), 
                  confirmed_case = mean(confirmed_case), 
                  dose = mean(dose)) %>%
        mutate(run_id = "mean_district")
      target_table_district <- plyr::rbind.fill(target_table_total_district_across_year, target_table_mean_district_across_year)
      # target_table <- plyr::rbind.fill(target_table, target_table_total_district_across_year, target_table_mean_district_across_year)
      # rm(target_table_total_district_across_year, target_table_mean_district_across_year)
    }

  }
  rm(target_table)
  

  # First the by year figure -- true case (non-cumulative and country-level mean)
  if(grepl("non", cumulative_type) & grepl("true", case_type)){
    plt_table <- target_table_country %>% 
      filter(threshold == chosen_threshold) %>%
      filter(general_scenario == "campaign_default" | surveillance_scenario == no_vaccination_surveillance_scenario) 
    coeff <- max(plt_table$dose) / max(plt_table$true_case)
    plt <- plt_table %>% 
      ggplot(aes(x=year, y=true_case, group=run_id)) +
      geom_line(data = plt_table[grepl("total_country", plt_table$run_id), ], 
                aes(x=year, y=true_case, group=run_id), 
                color = "lightblue", size = 0.05) +
      geom_line(data = plt_table[plt_table$run_id == "mean_country", ], 
                aes(x=year, y=true_case, group=run_id), 
                color = "red", size = 0.3) +
      geom_bar( data = plt_table[plt_table$run_id == "mean_country", ], 
                aes(x=year, y=dose / coeff), 
                stat="identity", size=.1, fill="#69b3a2", color="black", alpha=.4) + 
      scale_y_continuous( name = "True Cases",
                          sec.axis = sec_axis(trans=~.*coeff, name="Doses Administered")) + 
      facet_grid( ISO + admin_level ~ general_scenario + surveillance_scenario, scales = "free", space = "free") + 
      theme_minimal()
    # pdf("/home/kaiyuezou/VIMC_Model/surveillance_project/gavi_vimc_cholera/diagnostics/surveillance_project/test1.pdf")
    # plt
    # dev.off()
  }else if(grepl("non", cumulative_type) & grepl("clinical", case_type)){
    plt_table <- target_table_country %>% 
      filter(threshold == chosen_threshold) %>%
      filter(general_scenario == "campaign_default" | surveillance_scenario == no_vaccination_surveillance_scenario) 
    coeff <- max(plt_table$dose) / max(plt_table$clinical_case)
    plt <- plt_table %>% 
      ggplot(aes(x=year, y=clinical_case, group=run_id)) +
      geom_line(data = plt_table[grepl("total_country", plt_table$run_id), ], 
                aes(x=year, y=clinical_case, group=run_id), 
                color = "lightblue", size = 0.05) +
      geom_line(data = plt_table[plt_table$run_id == "mean_country", ], 
                aes(x=year, y=clinical_case, group=run_id), 
                color = "red", size = 0.3) +
      geom_bar( data = plt_table[plt_table$run_id == "mean_country", ], 
                aes(x=year, y=dose / coeff), 
                stat="identity", size=.1, fill="#69b3a2", color="black", alpha=.4) + 
      scale_y_continuous( name = "Clinical Cases",
                          sec.axis = sec_axis(trans=~.*coeff, name="Doses Administered")) + 
      facet_grid( ISO + admin_level ~ general_scenario + surveillance_scenario, scales = "free", space = "free") + 
      theme_minimal()
    # pdf("/home/kaiyuezou/VIMC_Model/surveillance_project/gavi_vimc_cholera/diagnostics/surveillance_project/test2.pdf")
    # plt
    # dev.off()
  }else if(grepl("non", cumulative_type) & grepl("confirmed", case_type)){
    plt_table <- target_table_country %>% 
      filter(threshold == chosen_threshold) %>%
      filter(general_scenario == "campaign_default" | surveillance_scenario == no_vaccination_surveillance_scenario) 
    coeff <- max(plt_table$dose) / max(plt_table$confirmed_case)
    plt <- plt_table %>% 
      ggplot(aes(x=year, y=confirmed_case, group=run_id)) +
      geom_line(data = plt_table[grepl("total_country", plt_table$run_id), ], 
                aes(x=year, y=confirmed_case, group=run_id), 
                color = "lightblue", size = 0.05) +
      geom_line(data = plt_table[plt_table$run_id == "mean_country", ], 
                aes(x=year, y=confirmed_case, group=run_id), 
                color = "red", size = 0.3) +
      geom_bar( data = plt_table[plt_table$run_id == "mean_country", ], 
                aes(x=year, y=dose / coeff), 
                stat="identity", size=.1, fill="#69b3a2", color="black", alpha=.4) + 
      scale_y_continuous( name = "Confirmed Cases",
                          sec.axis = sec_axis(trans=~.*coeff, name="Doses Administered")) + 
      facet_grid( ISO + admin_level ~ general_scenario + surveillance_scenario, scales = "free", space = "free") + 
      theme_minimal()
    # pdf("/home/kaiyuezou/VIMC_Model/surveillance_project/gavi_vimc_cholera/diagnostics/surveillance_project/test2.pdf")
    # plt
    # dev.off()
  }else if(grepl("non", cumulative_type) & grepl("averted_tr", case_type)){
    plt_table <- target_table_country %>% 
      filter(threshold == chosen_threshold) 
    coeff <- max(plt_table$dose) / max(plt_table$averted_true_case)
    plt <- plt_table %>% 
      ggplot(aes(x=year, y=averted_true_case, group=run_id)) +
      geom_line(data = plt_table[grepl("total_country", plt_table$run_id), ], 
                aes(x=year, y=averted_true_case, group=run_id), 
                color = "lightblue", size = 0.05) +
      geom_line(data = plt_table[plt_table$run_id == "mean_country", ], 
                aes(x=year, y=averted_true_case, group=run_id), 
                color = "red", size = 0.3) +
      geom_bar( data = plt_table[plt_table$run_id == "mean_country", ], 
                aes(x=year, y=dose / coeff), 
                stat="identity", size=.1, fill="#69b3a2", color="black", alpha=.4) + 
      scale_y_continuous( name = "Averted True Cases",
                          sec.axis = sec_axis(trans=~.*coeff, name="Doses Administered")) + 
      facet_grid( ISO + admin_level ~ surveillance_scenario, scales = "free", space = "free") + 
      theme_minimal()
    # pdf("/home/kaiyuezou/VIMC_Model/surveillance_project/gavi_vimc_cholera/diagnostics/surveillance_project/test3.pdf")
    # plt
    # dev.off()
  }else if(grepl("non", cumulative_type) & grepl("averted_cl", case_type)){
    plt_table <- target_table_country %>% 
      filter(threshold == chosen_threshold)
    coeff <- max(plt_table$dose) / max(plt_table$averted_clinical_case)
    plt <- plt_table %>% 
      ggplot(aes(x=year, y=averted_clinical_case, group=run_id)) +
      geom_line(data = plt_table[grepl("total_country", plt_table$run_id), ], 
                aes(x=year, y=averted_clinical_case, group=run_id), 
                color = "lightblue", size = 0.05) +
      geom_line(data = plt_table[plt_table$run_id == "mean_country", ], 
                aes(x=year, y=averted_clinical_case, group=run_id), 
                color = "red", size = 0.3) +
      geom_bar( data = plt_table[plt_table$run_id == "mean_country", ], 
                aes(x=year, y=dose / coeff), 
                stat="identity", size=.1, fill="#69b3a2", color="black", alpha=.4) + 
      scale_y_continuous( name = "Averted Clinical Cases",
                          sec.axis = sec_axis(trans=~.*coeff, name="Doses Administered")) + 
      facet_grid( ISO + admin_level ~ surveillance_scenario, scales = "free", space = "free") + 
      theme_minimal()
    # pdf("/home/kaiyuezou/VIMC_Model/surveillance_project/gavi_vimc_cholera/diagnostics/surveillance_project/test4.pdf")
    # plt
    # dev.off()
  }else if(grepl("non", cumulative_type) & grepl("averted_cf", case_type)){
    plt_table <- target_table_country %>% 
      filter(threshold == chosen_threshold)
    coeff <- max(plt_table$dose) / max(plt_table$averted_confirmed_case)
    plt <- plt_table %>% 
      ggplot(aes(x=year, y=averted_confirmed_case, group=run_id)) +
      geom_line(data = plt_table[grepl("total_country", plt_table$run_id), ], 
                aes(x=year, y=averted_confirmed_case, group=run_id), 
                color = "lightblue", size = 0.05) +
      geom_line(data = plt_table[plt_table$run_id == "mean_country", ], 
                aes(x=year, y=averted_confirmed_case, group=run_id), 
                color = "red", size = 0.3) +
      geom_bar( data = plt_table[plt_table$run_id == "mean_country", ], 
                aes(x=year, y=dose / coeff), 
                stat="identity", size=.1, fill="#69b3a2", color="black", alpha=.4) + 
      scale_y_continuous( name = "Averted Confirmed Cases",
                          sec.axis = sec_axis(trans=~.*coeff, name="Doses Administered")) + 
      facet_grid( ISO + admin_level ~ surveillance_scenario, scales = "free", space = "free") + 
      theme_minimal()
    # pdf("/home/kaiyuezou/VIMC_Model/surveillance_project/gavi_vimc_cholera/diagnostics/surveillance_project/test4.pdf")
    # plt
    # dev.off()
  }


  # Second: country-level cumulative across years 
  if(grepl("country", cumulative_type) & grepl("true", case_type)){
    plt_table <- target_table_country %>% 
      filter(threshold == chosen_threshold) %>%
      filter(general_scenario == "campaign_default" | surveillance_scenario == no_vaccination_surveillance_scenario) %>% 
      group_by(ISO, admin_level, incid_trend, outbk_trend, general_scenario, surveillance_scenario, run_id) %>%
      mutate(true_case = cumsum(true_case))
    plt <- plt_table %>% 
      ggplot(aes(x=year, y=true_case, group=run_id)) +
      geom_line(data = plt_table[plt_table$run_id == "mean_country", ], 
                aes(x=year, y=true_case, group=run_id), 
                color = "red", size = 0.3) +
      geom_line(data = plt_table[grepl("total_country", plt_table$run_id), ], 
                aes(x=year, y=true_case, group=run_id), 
                color = "lightblue", size = 0.05) +
      facet_grid( ISO + admin_level ~ general_scenario + surveillance_scenario, scales = "free", space = "free") + 
      theme_minimal()
    # pdf("/home/kaiyuezou/VIMC_Model/surveillance_project/gavi_vimc_cholera/diagnostics/surveillance_project/test5.pdf")
    # plt
    # dev.off()
  }else if(grepl("country", cumulative_type) & grepl("clinical", case_type)){
    plt_table <- target_table_country %>% 
      filter(threshold == chosen_threshold) %>%
      filter(general_scenario == "campaign_default" | surveillance_scenario == no_vaccination_surveillance_scenario) %>% 
      group_by(ISO, admin_level, incid_trend, outbk_trend, general_scenario, surveillance_scenario, run_id) %>%
      mutate(clinical_case = cumsum(clinical_case))
    plt <- plt_table %>% 
      ggplot(aes(x=year, y=clinical_case, group=run_id)) +
      geom_line(data = plt_table[plt_table$run_id == "mean_country", ], 
                aes(x=year, y=clinical_case, group=run_id), 
                color = "red", size = 0.3) +
      geom_line(data = plt_table[grepl("total_country", plt_table$run_id), ], 
                aes(x=year, y=clinical_case, group=run_id), 
                color = "lightblue", size = 0.05) +
      facet_grid( ISO + admin_level ~ general_scenario + surveillance_scenario, scales = "free", space = "free") + 
      theme_minimal()
    # pdf("/home/kaiyuezou/VIMC_Model/surveillance_project/gavi_vimc_cholera/diagnostics/surveillance_project/test6.pdf")
    # plt
    # dev.off()
  }else if(grepl("country", cumulative_type) & grepl("confirmed", case_type)){
    plt_table <- target_table_country %>% 
      filter(threshold == chosen_threshold) %>%
      filter(general_scenario == "campaign_default" | surveillance_scenario == no_vaccination_surveillance_scenario) %>% 
      group_by(ISO, admin_level, incid_trend, outbk_trend, general_scenario, surveillance_scenario, run_id) %>%
      mutate(confirmed_case = cumsum(confirmed_case))
    plt <- plt_table %>% 
      ggplot(aes(x=year, y=confirmed_case, group=run_id)) +
      geom_line(data = plt_table[plt_table$run_id == "mean_country", ], 
                aes(x=year, y=confirmed_case, group=run_id), 
                color = "red", size = 0.3) +
      geom_line(data = plt_table[grepl("total_country", plt_table$run_id), ], 
                aes(x=year, y=confirmed_case, group=run_id), 
                color = "lightblue", size = 0.05) +
      facet_grid( ISO + admin_level ~ general_scenario + surveillance_scenario, scales = "free", space = "free") + 
      theme_minimal()
    # pdf("/home/kaiyuezou/VIMC_Model/surveillance_project/gavi_vimc_cholera/diagnostics/surveillance_project/test6.pdf")
    # plt
    # dev.off()
  }else if(grepl("country", cumulative_type) & grepl("averted_tr", case_type)){
    plt_table <- target_table_country %>% 
      filter(threshold == chosen_threshold) %>% 
      group_by(ISO, admin_level, incid_trend, outbk_trend, surveillance_scenario, run_id) %>%
      mutate(averted_true_case = cumsum(averted_true_case))
    plt <- plt_table %>% 
      ggplot(aes(x=year, y=averted_true_case, group=run_id)) +
      geom_line(data = plt_table[plt_table$run_id == "mean_country", ], 
                aes(x=year, y=averted_true_case, group=run_id), 
                color = "red", size = 0.3) +
      geom_line(data = plt_table[grepl("total_country", plt_table$run_id), ], 
                aes(x=year, y=averted_true_case, group=run_id), 
                color = "lightblue", size = 0.05) +
      facet_grid( ISO + admin_level ~ surveillance_scenario, scales = "free", space = "free") + 
      theme_minimal()
    # pdf("/home/kaiyuezou/VIMC_Model/surveillance_project/gavi_vimc_cholera/diagnostics/surveillance_project/test7.pdf")
    # plt
    # dev.off()
  }else if(grepl("country", cumulative_type) & grepl("averted_cl", case_type)){
    plt_table <- target_table_country %>% 
      filter(threshold == chosen_threshold) %>% 
      group_by(ISO, admin_level, incid_trend, outbk_trend, surveillance_scenario, run_id) %>%
      mutate(averted_clinical_case = cumsum(averted_clinical_case))
    plt <- plt_table %>% 
      ggplot(aes(x=year, y=averted_clinical_case, group=run_id)) +
      geom_line(data = plt_table[plt_table$run_id == "mean_country", ], 
                aes(x=year, y=averted_clinical_case, group=run_id), 
                color = "red", size = 0.3) +
      geom_line(data = plt_table[grepl("total_country", plt_table$run_id), ], 
                aes(x=year, y=averted_clinical_case, group=run_id), 
                color = "lightblue", size = 0.05) +
      facet_grid( ISO + admin_level ~ surveillance_scenario, scales = "free", space = "free") + 
      theme_minimal()
    # pdf("/home/kaiyuezou/VIMC_Model/surveillance_project/gavi_vimc_cholera/diagnostics/surveillance_project/test8.pdf")
    # plt
    # dev.off()
  }else if(grepl("country", cumulative_type) & grepl("averted_cf", case_type)){
    plt_table <- target_table_country %>% 
      filter(threshold == chosen_threshold) %>% 
      group_by(ISO, admin_level, incid_trend, outbk_trend, surveillance_scenario, run_id) %>%
      mutate(averted_confirmed_case = cumsum(averted_confirmed_case))
    plt <- plt_table %>% 
      ggplot(aes(x=year, y=averted_confirmed_case, group=run_id)) +
      geom_line(data = plt_table[plt_table$run_id == "mean_country", ], 
                aes(x=year, y=averted_confirmed_case, group=run_id), 
                color = "red", size = 0.3) +
      geom_line(data = plt_table[grepl("total_country", plt_table$run_id), ], 
                aes(x=year, y=averted_confirmed_case, group=run_id), 
                color = "lightblue", size = 0.05) +
      facet_grid( ISO + admin_level ~ surveillance_scenario, scales = "free", space = "free") + 
      theme_minimal()
    # pdf("/home/kaiyuezou/VIMC_Model/surveillance_project/gavi_vimc_cholera/diagnostics/surveillance_project/test8.pdf")
    # plt
    # dev.off()
  }


  # The third: cumulative histogram by districts 
  if(grepl("district", cumulative_type) & grepl("true", case_type)){
    plt_table <- target_table_mean_district_across_year %>% 
      filter(threshold == chosen_threshold) %>%
      filter(general_scenario == "campaign_default" | surveillance_scenario == no_vaccination_surveillance_scenario)
    plt <- plt_table %>% 
      ggplot() +
      geom_bar(data = plt_table, 
              aes(x=NAME_1, y=true_case), 
              stat = "identity", color = "black", size = 0.2, fill = "orange") +
      facet_grid( ISO + admin_level ~ general_scenario + surveillance_scenario, scales = "free", space = "free") +
      theme_minimal() + 
      coord_flip()
    # pdf("/home/kaiyuezou/VIMC_Model/surveillance_project/gavi_vimc_cholera/diagnostics/surveillance_project/test9.pdf")
    # plt
    # dev.off()
  }else if(grepl("district", cumulative_type) & grepl("clinical", case_type)){
    plt_table <- target_table_mean_district_across_year %>% 
      filter(threshold == chosen_threshold) %>%
      filter(general_scenario == "campaign_default" | surveillance_scenario == no_vaccination_surveillance_scenario)
    plt <- plt_table %>% 
      ggplot() +
      geom_bar(data = plt_table, 
              aes(x=NAME_1, y=clinical_case), 
              stat = "identity", color = "black", size = 0.2, fill = "orange") +
      facet_grid( ISO + admin_level ~ general_scenario + surveillance_scenario, scales = "free", space = "free") +
      theme_minimal() + 
      coord_flip()
    # pdf("/home/kaiyuezou/VIMC_Model/surveillance_project/gavi_vimc_cholera/diagnostics/surveillance_project/test10.pdf")
    # plt
    # dev.off()
  }else if(grepl("district", cumulative_type) & grepl("confirmed", case_type)){
    plt_table <- target_table_mean_district_across_year %>% 
      filter(threshold == chosen_threshold) %>%
      filter(general_scenario == "campaign_default" | surveillance_scenario == no_vaccination_surveillance_scenario)
    plt <- plt_table %>% 
      ggplot() +
      geom_bar(data = plt_table, 
              aes(x=NAME_1, y=confirmed_case), 
              stat = "identity", color = "black", size = 0.2, fill = "orange") +
      facet_grid( ISO + admin_level ~ general_scenario + surveillance_scenario, scales = "free", space = "free") +
      theme_minimal() + 
      coord_flip()
    # pdf("/home/kaiyuezou/VIMC_Model/surveillance_project/gavi_vimc_cholera/diagnostics/surveillance_project/test10.pdf")
    # plt
    # dev.off()
  }else if(grepl("district", cumulative_type) & grepl("averted_tr", case_type)){
    plt_table <- target_table_mean_district_across_year %>% 
      filter(threshold == chosen_threshold)
    plt <- plt_table %>% 
      ggplot() +
      geom_bar(data = plt_table, 
              aes(x=NAME_1, y=averted_true_case), 
              stat = "identity", color = "black", size = 0.2, fill = "orange") +
      facet_grid( ISO + admin_level ~ surveillance_scenario, scales = "free", space = "free") +
      theme_minimal() + 
      coord_flip()
    # pdf("/home/kaiyuezou/VIMC_Model/surveillance_project/gavi_vimc_cholera/diagnostics/surveillance_project/test11.pdf")
    # plt
    # dev.off()
  }else if(grepl("district", cumulative_type) & grepl("averted_cl", case_type)){
    plt_table <- target_table_mean_district_across_year %>% 
      filter(threshold == chosen_threshold)
    plt <- plt_table %>% 
      ggplot() +
      geom_bar(data = plt_table, 
              aes(x=NAME_1, y=averted_clinical_case, fill=ISO), 
              stat = "identity", color = "black", size = 0.2, fill = "orange") +
      facet_grid( ISO + admin_level ~ surveillance_scenario, scales = "free", space = "free") +
      theme_minimal() + 
      coord_flip()
    # pdf("/home/kaiyuezou/VIMC_Model/surveillance_project/gavi_vimc_cholera/diagnostics/surveillance_project/test12.pdf")
    # plt
    # dev.off()
  }else if(grepl("district", cumulative_type) & grepl("averted_cf", case_type)){
    plt_table <- target_table_mean_district_across_year %>% 
      filter(threshold == chosen_threshold)
    plt <- plt_table %>% 
      ggplot() +
      geom_bar(data = plt_table, 
              aes(x=NAME_1, y=averted_confirmed_case, fill=ISO), 
              stat = "identity", color = "black", size = 0.2, fill = "orange") +
      facet_grid( ISO + admin_level ~ surveillance_scenario, scales = "free", space = "free") +
      theme_minimal() + 
      coord_flip()
    # pdf("/home/kaiyuezou/VIMC_Model/surveillance_project/gavi_vimc_cholera/diagnostics/surveillance_project/test12.pdf")
    # plt
    # dev.off()
  }


  # Return one single plot 
  return(plt)

}



#' @export
#' @name plot_efficiency 
#' @title plot_efficiency 
#' @description plot the efficacy and return the ggplot object
#' @param cache the cached environment
#' @param compare_type type of the comparison
#' @param threshold threshold to plot  
#' @param cumulative_type either by year (non_cumulative) or cumulated across years
#' @return cached file names 
plot_efficiency  <- function(cache, compare_type, threshold, cumulative_type){

  # Avoid using the same name for threshold
  chosen_threshold <- threshold

  # Restructure of the target table -- no NA's originally because it's all true case
  if(!"target_table_eff_side" %in% names(cache)){
    target_table <- cache$target_table %>%
      mutate(general_scenario = "both", 
            campaign_default_clinical_case = campaign_default_true_case / true_confirm_rate, 
            no_vaccination_clinical_case = no_vaccination_true_case / true_confirm_rate) %>% 
      mutate(averted_true_case = no_vaccination_true_case - campaign_default_true_case, 
            averted_clinical_case = no_vaccination_clinical_case - campaign_default_clinical_case) %>%
      mutate(dose = actual_fvp * 2, 
            efficacy = averted_true_case / actual_fvp)
    target_table$efficacy[is.infinite(target_table$efficacy)] <- 0
    target_table$efficacy[is.na(target_table$efficacy)] <- 0
    cache$target_table_eff_side <- target_table
  }else{
    target_table <- cache$target_table_eff_side
  }
  
  if(grepl("subtract", compare_type)){
    if(!"target_table_eff_sub" %in% names(cache)){
      table_no <- target_table[target_table$surveillance_scenario == "no-estimate", ] %>% 
        rename(efficacy_no_estimate = efficacy) %>% 
        select(ISO, admin_level, NAME_1, NAME_2, incid_trend, outbk_trend, threshold, year, run_id, efficacy_no_estimate)
      target_table <- target_table[target_table$surveillance_scenario != "no-estimate", ]
      target_table <- dplyr::left_join(target_table, table_no, by = c("ISO", "admin_level", "NAME_1", "NAME_2", "incid_trend", "outbk_trend", "threshold", "year", "run_id")) %>% 
        mutate(efficacy = efficacy - efficacy_no_estimate) %>% 
        select(- efficacy_no_estimate)
      cache$target_table_eff_sub <- target_table
    }else{
      target_table <- cache$target_table_eff_sub
    }
  }


  # For side-by-side comparison and direct comparison 
  if(grepl("side", compare_type) & grepl("non", cumulative_type)){
    target_table_total_country <- target_table %>% 
      group_by(ISO, admin_level, incid_trend, outbk_trend, general_scenario, surveillance_scenario, threshold, year, run_id) %>% 
      summarise(dose = sum(dose), 
                actual_fvp = sum(actual_fvp), 
                averted_true_case = sum(averted_true_case)) %>%
      mutate(run_id = paste0("total_country", run_id), efficacy = averted_true_case / actual_fvp)
    target_table_total_country$efficacy[is.infinite(target_table_total_country$efficacy)] <- 0
    target_table_total_country$efficacy[is.na(target_table_total_country$efficacy)] <- 0
    
  }else if(grepl("side", compare_type) & grepl("year", cumulative_type)){
    target_table_total_country <- target_table %>% 
      group_by(ISO, admin_level, incid_trend, outbk_trend, general_scenario, surveillance_scenario, threshold, run_id) %>% 
      summarise(dose = sum(dose), 
                actual_fvp = sum(actual_fvp), 
                averted_true_case = sum(averted_true_case)) %>%
      mutate(run_id = paste0("total_country", run_id), efficacy = averted_true_case / actual_fvp)
    target_table_total_country$efficacy[is.infinite(target_table_total_country$efficacy)] <- 0
    target_table_total_country$efficacy[is.na(target_table_total_country$efficacy)] <- 0
    
  }else if(grepl("side", compare_type) & grepl("district", cumulative_type)){
    target_table_total_country <- target_table %>% 
      group_by(ISO, admin_level, NAME_1, NAME_2, incid_trend, outbk_trend, general_scenario, surveillance_scenario, threshold, run_id) %>% 
      summarise(dose = sum(dose), 
                actual_fvp = sum(actual_fvp), 
                averted_true_case = sum(averted_true_case)) %>%
      mutate(efficacy = averted_true_case / actual_fvp)
    target_table_total_country$efficacy[is.infinite(target_table_total_country$efficacy)] <- 0
    target_table_total_country$efficacy[is.na(target_table_total_country$efficacy)] <- 0
    
  }else if(grepl("subtract", compare_type) & grepl("non", cumulative_type)){
    target_table_total_country <- target_table %>% 
      group_by(ISO, admin_level, incid_trend, outbk_trend, general_scenario, surveillance_scenario, threshold, year, run_id) %>% 
      summarise(efficacy = weighted.mean(efficacy, actual_fvp, na.rm = T)) 
    target_table_total_country$efficacy[is.infinite(target_table_total_country$efficacy)] <- 0
    target_table_total_country$efficacy[is.na(target_table_total_country$efficacy)] <- 0
  
  }else if(grepl("subtract", compare_type) & grepl("year", cumulative_type)){
    target_table_total_country <- target_table %>% 
      group_by(ISO, admin_level, incid_trend, outbk_trend, general_scenario, surveillance_scenario, threshold, run_id) %>% 
      summarise(efficacy = weighted.mean(efficacy, actual_fvp, na.rm = T)) 
    target_table_total_country$efficacy[is.infinite(target_table_total_country$efficacy)] <- 0
    target_table_total_country$efficacy[is.na(target_table_total_country$efficacy)] <- 0
  
  }else if(grepl("subtract", compare_type) & grepl("district", cumulative_type)){
    target_table_total_country <- target_table %>% 
      group_by(ISO, admin_level, NAME_1, NAME_2, incid_trend, outbk_trend, general_scenario, surveillance_scenario, threshold, run_id) %>% 
      summarise(efficacy = weighted.mean(efficacy, actual_fvp, na.rm = T)) 
    target_table_total_country$efficacy[is.infinite(target_table_total_country$efficacy)] <- 0
    target_table_total_country$efficacy[is.na(target_table_total_country$efficacy)] <- 0
    
  }


  # Plot
  if(grepl("side", compare_type) & grepl("non", cumulative_type)){
    plt_table <- target_table_total_country %>% 
      filter(threshold == chosen_threshold) %>% 
      mutate(year = as.factor(year))
    plt <- plt_table %>% 
      ggplot(aes(x=year, y=efficacy * 1000, fill=year)) +
      geom_boxplot() + 
      ylab("Averted true cases per 1000 fvps") + 
      facet_grid( ISO + admin_level ~ surveillance_scenario, scales = "free", space = "free") + 
      theme_minimal() + 
      coord_flip()
    # pdf("/home/kaiyuezou/VIMC_Model/surveillance_project/gavi_vimc_cholera/diagnostics/surveillance_project/test13.pdf")
    # plt
    # dev.off()
  }else if(grepl("side", compare_type) & grepl("year", cumulative_type)){
    plt_table <- target_table_total_country %>% 
      filter(threshold == chosen_threshold)
    plt <- plt_table %>% 
      ggplot(aes(y=efficacy * 1000, fill=ISO)) +
      geom_boxplot() + 
      ylab("Averted true cases per 1000 fvps") + 
      facet_grid( ISO + admin_level ~ surveillance_scenario, scales = "free", space = "free") + 
      theme_minimal()
    # pdf("/home/kaiyuezou/VIMC_Model/surveillance_project/gavi_vimc_cholera/diagnostics/surveillance_project/test14.pdf")
    # plt
    # dev.off()
  }else if(grepl("side", compare_type) & grepl("district", cumulative_type)){
    plt_table <- target_table_total_country %>% 
      filter(threshold == chosen_threshold)
    plt <- plt_table %>% 
      ggplot(aes(x=NAME_1, y=efficacy * 1000, fill=ISO)) +
      geom_boxplot() + 
      ylab("Averted true cases per 1000 fvps") + 
      facet_grid( ISO + admin_level ~ surveillance_scenario, scales = "free", space = "free") + 
      theme_minimal() + 
      coord_flip()
    # pdf("/home/kaiyuezou/VIMC_Model/surveillance_project/gavi_vimc_cholera/diagnostics/surveillance_project/test17.pdf", height = 20, width = 15)
    # plt
    # dev.off()
  }else if(grepl("subtract", compare_type) & grepl("non", cumulative_type)){
    plt_table <- target_table_total_country %>% 
      filter(threshold == chosen_threshold) %>% 
      mutate(year = as.factor(year))
    plt <- plt_table %>% 
      ggplot(aes(x=year, y=efficacy * 1000, fill=year)) +
      geom_boxplot() + 
      ylab("Averted true cases per 1000 fvps") + 
      facet_grid( ISO + admin_level ~ surveillance_scenario, scales = "free", space = "free") + 
      theme_minimal() + 
      coord_flip()
    # pdf("/home/kaiyuezou/VIMC_Model/surveillance_project/gavi_vimc_cholera/diagnostics/surveillance_project/test15.pdf")
    # plt
    # dev.off()
  }else if(grepl("subtract", compare_type) & grepl("year", cumulative_type)){
    plt_table <- target_table_total_country %>% 
      filter(threshold == chosen_threshold)
    plt <- plt_table %>% 
      ggplot(aes(y=efficacy * 1000, fill=ISO)) +
      geom_boxplot() + 
      ylab("Averted true cases per 1000 fvps") + 
      facet_grid( ISO + admin_level ~ surveillance_scenario, scales = "free", space = "free") + 
      theme_minimal()
    # pdf("/home/kaiyuezou/VIMC_Model/surveillance_project/gavi_vimc_cholera/diagnostics/surveillance_project/test16.pdf")
    # plt
    # dev.off()
  }else if(grepl("subtract", compare_type) & grepl("district", cumulative_type)){
    plt_table <- target_table_total_country %>% 
      filter(threshold == chosen_threshold)
    plt <- plt_table %>% 
      ggplot(aes(x=NAME_1, y=efficacy * 1000, fill=ISO)) +
      geom_boxplot() + 
      ylab("Averted true cases per 1000 fvps") + 
      facet_grid( ISO + admin_level ~ surveillance_scenario, scales = "free", space = "free") + 
      theme_minimal() + 
      coord_flip()
    # pdf("/home/kaiyuezou/VIMC_Model/surveillance_project/gavi_vimc_cholera/diagnostics/surveillance_project/test18.pdf", height = 20, width = 15)
    # plt
    # dev.off()
  }


  # Return one single plot 
  return(plt)

}



#' @export
#' @name plot_time_table
#' @title plot_time_table
#' @description plot the time table
#' @param cache the cached environment
#' @return the kable object 
plot_time_table <- function(cache){
  time_table <- cache$time_table %>%
    dplyr::select(country, general_scenario, surveillance_scenario, threshold, year, 
                  load_baseline_incidence, update_targets_list, update_save_vac_raster, 
                  update_save_sus_raster, create_expectedCases, add_new_row_to_target_list) 
  time_table %>% 
    # dplyr::mutate_if(is.numeric, function(x) {format(x , big.mark=",")}) %>%
    kableExtra::kable(col.names = names(time_table)) %>%
    kableExtra::kable_styling(bootstrap_options = c("striped"), fixed_thead = T) %>%
    kableExtra::kable_paper(full_width = F) %>%
    kableExtra::row_spec(0, bold = T)
}



#' @name get_inc_tp_doses_helper
#' @title get_inc_tp_doses_helper
#' @description Function to assist calculation of get_inc_tp_doses(), read in only one data frame a time
#' @param rc a data frame of one scenario of one country
#' @param sum_level takes either "country" or "admin". If takes "admin", returns a table by admin units; if takes "country", returns a table for the whole country aggregating all admin units.
#' @param vax_cov vaccine coverage, default value is 0.68
#' @param surveillance_scenario scale of surveillance scenario, which will be used to calculate observed incidence 
#' @return a table with incidence, target population and number of doses administered, by year and cumulatively
#' @export
#' @include
get_inc_tp_doses_helper <- function(rc,
                                    sum_level,
                                    vax_cov = 0.68,
                                    surveillance_scenario){
  
  country <- rc$ISO[1]
  vax_years <- unique((rc %>% filter(!is.na(confirmation_lens)))$year)
  
  # fill in missing observed incidence column
  if(surveillance_scenario == "district-estimate"){
    rc <- rc %>% mutate(observed_incidence_rate = true_incidence_rate)
  }
  if(surveillance_scenario == "global-estimate"){
    global_comfirmation_rate <- unique(rc$confirmation_rate)[!is.na(unique(rc$confirmation_rate))]
    rc <- rc %>% mutate(observed_incidence_rate = true_incidence_rate / global_comfirmation_rate)
  }
  if(surveillance_scenario == "no-estimate"){
    rc <- rc %>% mutate(observed_incidence_rate = true_incidence_rate / confirmation_rate)
  }
  

  ### if want to summarize all admins together 
  if(sum_level == "country"){
    
    # by year
    df_tp_byyear <- rc %>% 
      filter(is_target == 1) %>%
      group_by(year, run_id) %>%
      summarize(target_pop = sum(pop_model)) %>%
      group_by(year) %>%
      summarize(tp_lb = quantile(target_pop, 0.025, na.rm = TRUE), # lower 95%CI
                tp_median = median(target_pop, , na.rm = TRUE),
                tp_ub = quantile(target_pop, 0.975, na.rm = TRUE)) %>% # upper 95%CI
      mutate(ISO = country) %>%
      dplyr::select(ISO, year, tp_lb, tp_median, tp_ub)
    
    
    # cumulatively
    df_tp_cumu <- rc %>% 
      filter(is_target == 1) %>%
      group_by(year, run_id) %>%
      summarize(target_pop = sum(pop_model)) %>%
      group_by(run_id) %>%
      mutate(target_pop_cumu = cumsum(target_pop)) %>%
      dplyr::select(-target_pop) %>%
      group_by(year) %>%
      summarize(tp_cumu_lb = quantile(target_pop_cumu, 0.025, na.rm = TRUE), # lower 95%CI
                tp_cumu_median = median(target_pop_cumu, na.rm = TRUE),
                tp_cumu_ub = quantile(target_pop_cumu, 0.975, na.rm = TRUE)) %>% # upper 95%CI
      mutate(ISO = country) %>%
      dplyr::select(ISO, year, tp_cumu_lb, tp_cumu_median, tp_cumu_ub)
    
    # join by year and cumulative tables 
    df_tp <- left_join(df_tp_byyear, df_tp_cumu, by = c("ISO", "year"))
    
    # fill in no targeting years with 0s
    df_add <- data.frame("ISO" = country, "year" = vax_years[vax_years %in% df_tp$year == FALSE], 
                         "tp_lb" = 0, "tp_median" = 0, "tp_ub" = 0,
                         "tp_cumu_lb" = 0, "tp_cumu_median" = 0, "tp_cumu_ub" = 0) %>% as.data.frame()
    df_tp <- rbind(df_tp, df_add) %>% arrange(year)
    
    
    # add incidence columns
    df_case <- rc %>% 
      mutate(observed_case = pop_model * observed_incidence_rate) %>%
      group_by(year, run_id) %>%
      summarize(observed_case = sum(observed_case),
                true_case = sum(true_case)) %>%
      group_by(year) %>%
      summarize(observed_case_lb = quantile(observed_case, 0.025, na.rm = TRUE),
                observed_case_median = median(observed_case, na.rm = TRUE),
                observed_case_ub = quantile(observed_case, 0.975, na.rm = TRUE),
                true_case_lb = quantile(true_case, 0.025, na.rm = TRUE),
                true_case_median = median(true_case, na.rm = TRUE),
                true_case_ub = quantile(true_case, 0.975, na.rm = TRUE))
    
    pop_by_year <- rc %>%
      group_by(year, run_id) %>%
      summarize(pop=sum(pop_model)) %>%
      filter(run_id == 1) %>%
      dplyr::select(-run_id)
    
    df_inc <- df_case %>%
      left_join(pop_by_year, by = "year") %>%
      mutate(true_inc_rate_lb = true_case_lb / pop,
             true_inc_rate_median = true_case_median / pop,
             true_inc_rate_ub = true_case_ub / pop,
             obs_inc_rate_lb = observed_case_lb / pop,
             obs_inc_rate_median = observed_case_median / pop,
             obs_inc_rate_ub = observed_case_ub / pop) %>%
      dplyr::select(-c(true_case_lb, true_case_median, true_case_ub,
                       observed_case_lb, observed_case_median, observed_case_ub,
                       pop))
    
    # join df_tp with df_inc 
    df_result <- df_inc %>%
      left_join(df_tp, by = "year") %>%
      mutate(ISO = country) %>%
      dplyr::select(ISO, c(1:7, 9:14))
    
  }
  
  ### if want to get data for each admin
  if(sum_level == "admin"){
    
    if("NAME_2" %in% names(rc)){
      rc <- rc %>% dplyr::rename(admin_name = NAME_2)
      admin_level <- "admin2"
    }else{
      rc <- rc %>% dplyr::rename(admin_name = NAME_1)
      admin_level <- "admin1"
    }
    
    # target pop by year
    df_tp_byyear <- rc %>% 
      filter(is_target == 1) %>%
      group_by(admin_name, year) %>% 
      summarize(tp_lb = quantile(pop_model, 0.025, na.rm = TRUE), # lower 95%CI
                tp_median = median(pop_model, na.rm = TRUE),
                tp_ub = quantile(pop_model, 0.975, na.rm = TRUE)) %>% # upper 95%CI
      mutate(ISO = country) %>%
      dplyr::select(ISO, admin_name, year, tp_lb, tp_median, tp_ub) 

  
    # target pop cumulatively 
    df_tp_cumu <- rc %>% 
      filter(is_target == 1) %>%
      group_by(admin_name, year, run_id) %>%
      summarize(target_pop = sum(pop_model)) %>%
      group_by(run_id) %>%
      mutate(target_pop_cumu = cumsum(target_pop)) %>%
      dplyr::select(-target_pop) %>%
      group_by(admin_name, year) %>%
      summarize(tp_cumu_lb = quantile(target_pop_cumu, 0.025, na.rm = TRUE), # lower 95%CI
                tp_cumu_median = median(target_pop_cumu, na.rm = TRUE),
                tp_cumu_ub = quantile(target_pop_cumu, 0.975, na.rm = TRUE)) %>% # upper 95%CI
      mutate(ISO = country) %>%
      dplyr::select(ISO, admin_name, year, tp_cumu_lb, tp_cumu_median, tp_cumu_ub)
      
    
    # join by year and cumulative tables 
    df_tp <- left_join(df_tp_byyear, df_tp_cumu, by = c("ISO", "admin_name", "year"))
    
    # add back all the admin units/years that are not targeted, value = 0
    df_temp <- data.frame("admin_name" = rep(unique(rc$admin_name), times = length(vax_years)),
                          "year" = rep(vax_years, each=length(unique(rc$admin_name))))
    df_tp <- left_join(df_temp, df_tp, by = c("admin_name", "year")) %>%
      mutate(ISO = country) %>%
      dplyr::select(ISO, admin_name, year, tp_lb, tp_median, tp_ub, tp_cumu_lb, tp_cumu_median, tp_cumu_ub)
    
    # replace NA with 0
    df_tp[is.na(df_tp)] <-0
    
    
    ## add incidence columns
    df_inc <- rc %>%
      group_by(admin_name, year) %>%
      summarize(true_inc_rate_lb = quantile(true_incidence_rate, 0.025, na.rm = TRUE), 
                true_inc_rate_median = median(true_incidence_rate, na.rm = TRUE),
                true_inc_rate_ub = quantile(true_incidence_rate, 0.975, na.rm = TRUE),
                obs_inc_rate_lb = quantile(observed_incidence_rate, 0.025, na.rm = TRUE), 
                obs_inc_rate_median = median(observed_incidence_rate, na.rm = TRUE),
                obs_inc_rate_ub = quantile(observed_incidence_rate, 0.975, na.rm = TRUE))
    
    # join df_inc and df_tp
    df_result <- df_inc %>%
      left_join(df_tp, by = c("admin_name", "year")) %>%
      mutate(ISO = country) %>%
      arrange(year, admin_name) %>%
      dplyr::select(ISO, admin_name, year, true_inc_rate_lb, true_inc_rate_median, true_inc_rate_ub, obs_inc_rate_lb, obs_inc_rate_median, obs_inc_rate_ub, tp_lb, tp_median, tp_ub, tp_cumu_lb, tp_cumu_median, tp_cumu_ub)
    colnames(df_result)[2] <- admin_level
    
  }
  
  # get number of doses administered
  df_result <- df_result %>%
    mutate(doses_lb = vax_cov * tp_lb * 2, doses_median = vax_cov * tp_median * 2, doses_ub = vax_cov * tp_ub * 2,
           doses_cumu_lb = vax_cov * tp_cumu_lb * 2, doses_cumu_median = vax_cov * tp_cumu_median * 2, doses_cumu_ub = vax_cov * tp_cumu_ub * 2)

  return(df_result)
}



#' @name get_inc_tp_doses
#' @title get_inc_tp_doses
#' @description function to get incidence, target population and number of doses administered under selected campaign scenario, surveillance scenario and threshold
#' @param rawoutpath path of output_raw in which all the model output files are stored
#' @param configpath path of configuration files
#' @param countries vector of all the countries included the diagnostic report
#' @param scenario takes value of either "campaign-default" or "no-vaccination"
#' @param threshold threshold chosen for the diagnostic report, takes one of 1/1000, 1/5000, or 1/10,000
#' @param surveillance_scenario scale of surveillance, takes "no-estimate", "global-estimate", or "district-estimate"
#' @return 
#' @export
#' @include
get_inc_tp_doses <- function(rawoutpath = "output_raw/202110gavi-3",
                             configpath = "configs/202110gavi-3",
                             countries,
                             scenario,
                             surveillance_scenario,
                             threshold,
                             admin_level, 
                             ...){
    
    # first, check whether desired output files are all in the output_raw folder 
    outputs_availability <- check_outputs_availability(rawoutpath = rawoutpath ,configpath = configpath, countries = countries,
                               scenario = scenario, surveillance_scenario = surveillance_scenario, threshold = threshold, admin_level = admin_level)
    if(outputs_availability == FALSE){
        stop(paste("Stop calculating target population as not all the required output files are available in", rawoutpath, "folder."))
    }else{

        # set up file path
        output_path <- paste0(rawoutpath, "/", scenario)
        df_setting <- check_model_setting(configpath = configpath, countries = countries,
                                          scenario = scenario, surveillance_scenario = surveillance_scenario, 
                                          threshold = threshold)
        
        # loop through countries
        for(country in countries){

            # get filenames for both admin levels
            incid <- df_setting[which(df_setting$country == country), ]$incidence_trend 
            outbk <- df_setting[which(df_setting$country == country), ]$outbreak_multiplier 
            vax_cov <- df_setting[which(df_setting$country == country), ]$vac_coverage
            surveillance_scenario <- df_setting[which(df_setting$country == country), ]$surveillance_scenario
            prefix <- paste0("incid_", incid, "_outbk_", outbk, "_")
            if(admin_level %in% c("both", "admin1")){filename_admin1 <- paste0(output_path, "/", prefix, threshold, "_", surveillance_scenario, "_", country, "_rc_admin1.csv")}
            if(admin_level %in% c("both", "admin2")){filename_admin2 <- paste0(output_path, "/", prefix, threshold, "_", surveillance_scenario, "_", country, "_rc_admin2.csv")}

            # read in csv files
            if(admin_level %in% c("both", "admin1")){df_admin1 <- read.csv(filename_admin1)}
            if(admin_level %in% c("both", "admin2")){df_admin2 <- read.csv(filename_admin2)}

            # if this is the first country, make new tables of target pop and doses administered
            if(which(countries == country) == 1){
                # when doing admin1 level OCV targeting
                if(admin_level %in% c("both", "admin1")){
                  df_admin1_sumadmins <- get_inc_tp_doses_helper(rc = df_admin1, sum_level = "country", vax_cov = vax_cov, surveillance_scenario = surveillance_scenario) # sum all admins together
                  df_admin1_byadmin <- get_inc_tp_doses_helper(rc = df_admin1, sum_level = "admin", vax_cov = vax_cov, surveillance_scenario = surveillance_scenario) %>% 
                    dplyr::arrange(admin1, year) # by admin
                }
                # when doing admin2 level OCV targeting
                if(admin_level %in% c("both", "admin2")){
                  df_admin2_sumadmins <- get_inc_tp_doses_helper(rc = df_admin2, sum_level = "country", vax_cov = vax_cov, surveillance_scenario = surveillance_scenario) # sum all admins together
                  df_admin2_byadmin <- get_inc_tp_doses_helper(rc = df_admin2, sum_level = "admin", vax_cov = vax_cov, surveillance_scenario = surveillance_scenario) %>% 
                    dplyr::arrange(admin2, year) # by admin
                }
            }else{ # else, append to existing result tables
                if(admin_level %in% c("both", "admin1")){
                  df_admin1_sumadmins <- rbind(df_admin1_sumadmins, get_inc_tp_doses_helper(rc = df_admin1, sum_level = "country", vax_cov = vax_cov, surveillance_scenario = surveillance_scenario))
                  df_admin1_byadmin <- rbind(df_admin1_byadmin, get_inc_tp_doses_helper(rc = df_admin1, sum_level = "admin", vax_cov = vax_cov, surveillance_scenario = surveillance_scenario)) %>% 
                    dplyr::arrange(admin1, year)
                }
                if(admin_level %in% c("both", "admin2")){
                  df_admin2_sumadmins <- rbind(df_admin2_sumadmins, get_inc_tp_doses_helper(rc = df_admin2, sum_level = "country", vax_cov = vax_cov, surveillance_scenario = surveillance_scenario))
                  df_admin2_byadmin <- rbind(df_admin2_byadmin, get_inc_tp_doses_helper(rc = df_admin2, sum_level = "admin", vax_cov = vax_cov, surveillance_scenario = surveillance_scenario)) %>% 
                    dplyr::arrange(admin2, year)
                }
            }
        }
    }
    # make a list storing all 4 result tables
    if(!exists("df_admin1_sumadmins")){df_admin1_sumadmins <- NULL}
    if(!exists("df_admin1_byadmin")){df_admin1_byadmin <- NULL}
    if(!exists("df_admin2_sumadmins")){df_admin2_sumadmins <- NULL}
    if(!exists("df_admin2_byadmin")){df_admin2_byadmin <- NULL}
    
    list_inc_tp_doses <- list("df_admin1_sumadmins" = df_admin1_sumadmins,
                              "df_admin1_byadmin" = df_admin1_byadmin,
                              "df_admin2_sumadmins" = df_admin2_sumadmins,
                              "df_admin2_byadmin" = df_admin2_byadmin)

    return(list_inc_tp_doses)
}


#' @name omega_alpha_distribution
#' @title omega_alpha_distribution
#' @description plot the distribution of omega d
#' @param countries vector of all the countries included the diagnostic report
#' @param surveillance_scenario_chosen scale of surveillance, takes "no-estimate", "global-estimate", or "district-estimate"
#' @param admins admin levels to include
#' @param year_chosen 
#' @param threshold_chosen 
#' @param cache where the filenames and files are saved 
#' @return two distribution plots cached 
#' @export
#' @include
omega_alpha_distribution <- function( countries,
                                      surveillance_scenario_chosen = "global-estimate",   
                                      admins = "both", 
                                      year_chosen = 2022, 
                                      threshold_chosen = 0.001, 
                                      cache, 
                                      ...){   
  # Get the distribution table
  omega_dist <- cache$target_table %>%
    filter(ISO %in% countries, confirmation_lens == surveillance_scenario_chosen, year == year_chosen, threshold == threshold_chosen) %>%
    {if(admins != "both") filter(., admin_level == admins)} %>%
    select(ISO, admin_level, run_id, true_confirm_rate)
  alpha_dist <- cache$target_table %>%
    filter(ISO %in% countries, confirmation_lens == surveillance_scenario_chosen, year == year_chosen, threshold == threshold_chosen) %>%
    {if(admins != "both") filter(., admin_level == admins)} %>%
    group_by(ISO, admin_level, run_id) %>% 
    summarize(confirmation_rate = confirmation_rate[1]) %>% 
    select(ISO, admin_level, run_id, confirmation_rate)
  cache$alpha_table <- alpha_dist
  
  # Plot them
  plt_omega <- omega_dist %>% 
    ggplot(aes(x=true_confirm_rate)) +
    geom_histogram(fill="#69b3a2", color="#e9ecef", alpha=0.9) +
    facet_grid( ISO ~ admin_level ) + 
    theme_minimal() + 
    xlim(0, 1) + 
    labs(x = "omega_d")
  # pdf("/home/kaiyuezou/VIMC_Model/surveillance_project/gavi_vimc_cholera/test_omega.pdf")
  # plt_omega
  # dev.off()

  plt_alpha <- alpha_dist %>% 
    ggplot(aes(x=confirmation_rate)) +
    geom_histogram(fill="#69b3a2", color="#e9ecef", alpha=0.9) +
    facet_grid( ISO ~ admin_level ) + 
    theme_minimal() + 
    xlim(0, 1) + 
    labs(x = "alpha")
  # pdf("/home/kaiyuezou/VIMC_Model/surveillance_project/gavi_vimc_cholera/test_alpha.pdf")
  # plt_alpha
  # dev.off()

  # Cache them 
  cache$plt_omega <- plt_omega
  cache$plt_alpha <- plt_alpha

}



#' @name alpha_table
#' @title alpha_table
#' @description make a table of alpha
#' @param cache where the filenames and files are saved 
#' @return alpha table
#' @export
#' @include
alpha_table <- function(cache){   
  alpha_table <- cache$alpha_table %>% 
    mutate(run_id = paste0("run ", run_id)) %>% 
    tidyr::spread(run_id, confirmation_rate)

  alpha_table %>% 
    dplyr::mutate_if(is.numeric, function(x) {round(x, 3)}) %>%
    kableExtra::kable(col.names = names(alpha_table)) %>%
    kableExtra::kable_styling(bootstrap_options = c("striped"), fixed_thead = T) %>%
    kableExtra::kable_paper(full_width = F) %>%
    kableExtra::row_spec(0, bold = T)

}


#' @name make_eff_table
#' @title make_eff_table
#' @description make a table summarizing (total) target population and ocv efficiency (at country level)
#' @return a table summarizing target population and ocv efficiency for one country at one admin_level (either admin1 or admin2)
#' @export
#' @include
make_eff_table <- function(rc,
                           vax_cov = 0.68){
  
  country <- rc$ISO[1]

    df_tp_cumu <- rc %>% 
      filter(is_target == 1) %>%
      group_by(year, run_id, threshold, confirmation_lens) %>%
      summarize(target_pop = sum(pop_model)) %>%
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
      # calculate efficiency
      mutate(efficiency = true_ac_cumu / fvp_cumu * 1000) %>%
      group_by(threshold, confirmation_lens) %>%
      summarize(efficiency_lb = quantile(efficiency, 0.025, na.rm = TRUE), # lower 95%CI
                efficiency_median = median(efficiency, na.rm = TRUE),
                efficiency_ub = quantile(efficiency, 0.975, na.rm = TRUE)) %>% # upper 95%CI
      mutate(ISO = country) %>%
      dplyr::select(ISO, threshold, confirmation_lens, efficiency_lb, efficiency_median, efficiency_ub)
    
    
    # join two tables and round up numbers
    df_result <- df_tp_cumu %>%
      left_join(df_eff, by = c("ISO", "threshold", "confirmation_lens")) %>%
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
      select(ISO, threshold, confirmation_lens, tp_cumu, efficiency)
  
  return(df_result)
}


#' @name combine_eff_table
#' @title combine_eff_table
#' @description combine tables from multiple countries
#' @return a table summarizing target population and ocv efficiency for one or more countries (at country level)
#' @export
#' @include
combine_eff_table <- function(countries, 
                              rc,
                              admin_level){
    
    for(country in countries){
      
      if(admin_level %in% c("both", "admin1")){
        df_admin1 <- df_target %>% filter(admin_level == "admin1" & ISO == country)
        eff_table_admin1 <- make_eff_table(df_admin1) %>% mutate(admin_level = "admin1")
      }
      if(admin_level %in% c("both", "admin2")){
        df_admin2 <- df_target %>% filter(admin_level == "admin2" & ISO == country)
        eff_table_admin2 <- make_eff_table(df_admin2) %>% mutate(admin_level = "admin2")
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

