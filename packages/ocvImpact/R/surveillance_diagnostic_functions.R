###################
# Functions
###################

#' @export
#' @name get_filenames
#' @title get_filenames
#' @description plot the polygon with modeled cases
#' @param cache the cached environment
#' @param surveillance_project_directory surveillance project directory
#' @param pre_country countries specified in the parameters 
#' @param pre_vac_incid_thresholds incidence rate threshold specified in the parameters 
#' @param no_vaccination_surveillance_scenario the surveillance scenario specified for the no vaccination simulation 
#' @return cached file names 
get_filenames <- function(cache, surveillance_project_directory, pre_country, pre_vac_incid_thresholds, 
                          no_vaccination_surveillance_scenario){
  
  # Get the set_all_parameters.R to have the context first
  source(file.path(surveillance_project_directory, "scripts/set_all_parameters.R"))

  
  # Get the correct version of country list and threshold list
  country_list <- ifelse(pre_country %in% countries, pre_country, countries)  
  rm(pre_country)

  threshold_list <- ifelse(pre_vac_incid_thresholds %in% vac_incid_thresholds, pre_vac_incid_thresholds, vac_incid_thresholds)
  rm(pre_vac_incid_thresholds)

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
    if(sum(!file.exists(ecraster_fns)) > 0){stop(paste("Incomplete model outputs for expected cases rasters for country", country_list[i]))}

    ## expected cases table
    ectable_fns_camp <- apply(expand.grid(paste0("output_raw/", runname, "/campaign-default/", country_list[i]), 
                              "ec", paste0(admin_list, ".csv")), 1, paste, collapse = "_")
    ectable_fns_novac <- apply(expand.grid(paste0("output_raw/", runname, "/no-vaccination/", country_list[i]), 
                              "ec", paste0(admin_list, ".csv")), 1, paste, collapse = "_")                          
    ectable_fns <- c(ectable_fns_camp, ectable_fns_novac)
    ectable_fns <- paste0(surveillance_project_directory, "/", ectable_fns)
    if(sum(!file.exists(ectable_fns)) > 0){stop(paste("Incomplete model outputs for expected cases table for country", country_list[i]))}

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
    if(sum(!file.exists(truecomfirmrate_fns)) > 0){stop(paste("Incomplete model outputs for true confirmation rate raster for country", country_list[i]))}                        

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
#' @return cached file names 
plot_cases <- function(cache, case_type, threshold, cumulative_type){
  
  # Restructure of the target table 
  if(grepl("true", case_type) | grepl("clinical", case_type)){
    if(!"target_table_long" %in% names(cache)){
      target_table <- cache$target_table %>% 
        rename(campaign_default = campaign_default_true_case) %>%
        rename(no_vaccination = no_vaccination_true_case) 
      target_table <- target_table %>% 
        tidyr::gather("general_scenario", "true_case", match(c("campaign_default", "no_vaccination"), names(target_table))) %>%
        mutate(clinical_case = true_case / true_confirm_rate)
      cache$target_table_long <- target_table
    }else{
      target_table <- cache$target_table_long
    }
    
  }else{
    target_table <- cache$target_table %>%
      mutate(general_scenario = "both", 
            campaign_default_clinical_case = campaign_default_true_case / true_confirm_rate, 
            no_vaccination_clinical_case = no_vaccination_true_case / true_confirm_rate) %>% 
      mutate(averted_true_case = no_vaccination_true_case - campaign_default_true_case, 
            averted_clinical_case = no_vaccination_clinical_case - campaign_default_clinical_case)
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
                  averted_true_case = sum(averted_true_case), 
                  averted_clinical_case = sum(averted_clinical_case)) %>%
        mutate(run_id = paste0("total_country", run_id))
      target_table_mean_country <- target_table_total_country %>%
        group_by(ISO, admin_level, incid_trend, outbk_trend, general_scenario, surveillance_scenario, threshold, year) %>% 
        summarise(campaign_default_true_case = mean(campaign_default_true_case), 
                  no_vaccination_true_case = mean(no_vaccination_true_case), 
                  campaign_default_clinical_case = mean(campaign_default_clinical_case), 
                  no_vaccination_clinical_case = mean(no_vaccination_clinical_case), 
                  averted_true_case = mean(averted_true_case), 
                  averted_clinical_case = mean(averted_clinical_case)) %>%
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
                  averted_true_case = sum(averted_true_case), 
                  averted_clinical_case = sum(averted_clinical_case)) %>%
        mutate(run_id = "total_district")
      target_table_mean_district_across_year <- target_table_total_district_across_year %>%
        group_by(ISO, admin_level, NAME_1, NAME_2, incid_trend, outbk_trend, general_scenario, surveillance_scenario, threshold) %>% 
        summarise(campaign_default_true_case = mean(campaign_default_true_case), 
                  no_vaccination_true_case = mean(no_vaccination_true_case), 
                  campaign_default_clinical_case = mean(campaign_default_clinical_case), 
                  no_vaccination_clinical_case = mean(no_vaccination_clinical_case), 
                  averted_true_case = mean(averted_true_case), 
                  averted_clinical_case = mean(averted_clinical_case)) %>%
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
                  clinical_case = sum(clinical_case)) %>%
        mutate(run_id = paste0("total_country", run_id))
      target_table_mean_country <- target_table_total_country %>%
        group_by(ISO, admin_level, incid_trend, outbk_trend, general_scenario, surveillance_scenario, threshold, year) %>% 
        summarise(true_case = mean(true_case), 
                  clinical_case = mean(clinical_case)) %>%
        mutate(run_id = "mean_country")
      target_table_country <- plyr::rbind.fill(target_table_total_country, target_table_mean_country)
      # target_table <- plyr::rbind.fill(target_table, target_table_total_country, target_table_mean_country)
      # rm(target_table_total_country, target_table_mean_country)
    }else{
      target_table_total_district_across_year <- target_table %>%
        group_by(ISO, admin_level, NAME_1, NAME_2, incid_trend, outbk_trend, general_scenario, surveillance_scenario, threshold, run_id) %>% 
        summarise(true_case = sum(true_case), 
                  clinical_case = sum(clinical_case)) %>%
        mutate(run_id = "total_district")
      target_table_mean_district_across_year <- target_table_total_district_across_year %>%
        group_by(ISO, admin_level, NAME_1, NAME_2, incid_trend, outbk_trend, general_scenario, surveillance_scenario, threshold) %>% 
        summarise(true_case = mean(true_case), 
                  clinical_case = mean(clinical_case)) %>%
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
      filter(threshold == threshold) %>%
      filter(general_scenario == "campaign_default" | surveillance_scenario == "district-estimate") 
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
    # pdf("/home/kaiyuezou/VIMC_Model/surveillance_project/gavi_vimc_cholera/diagnostics/surveillance_project/test1.pdf")
    # plt
    # dev.off()
  }else if(grepl("non", cumulative_type) & grepl("clinical", case_type)){
    plt_table <- target_table_country %>% 
      filter(threshold == threshold) %>%
      filter(general_scenario == "campaign_default" | surveillance_scenario == "district-estimate") 
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
    # pdf("/home/kaiyuezou/VIMC_Model/surveillance_project/gavi_vimc_cholera/diagnostics/surveillance_project/test2.pdf")
    # plt
    # dev.off()
  }else if(grepl("non", cumulative_type) & grepl("averted_tr", case_type)){
    plt_table <- target_table_country %>% 
      filter(threshold == threshold) 
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
    # pdf("/home/kaiyuezou/VIMC_Model/surveillance_project/gavi_vimc_cholera/diagnostics/surveillance_project/test3.pdf")
    # plt
    # dev.off()
  }else if(grepl("non", cumulative_type) & grepl("averted_ob", case_type)){
    plt_table <- target_table_country %>% 
      filter(threshold == threshold)
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
    # pdf("/home/kaiyuezou/VIMC_Model/surveillance_project/gavi_vimc_cholera/diagnostics/surveillance_project/test4.pdf")
    # plt
    # dev.off()
  }


  # Second: country-level cumulative across years 
  if(grepl("country", cumulative_type) & grepl("true", case_type)){
    plt_table <- target_table_country %>% 
      filter(threshold == threshold) %>%
      filter(general_scenario == "campaign_default" | surveillance_scenario == "district-estimate") %>% 
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
      filter(threshold == threshold) %>%
      filter(general_scenario == "campaign_default" | surveillance_scenario == "district-estimate") %>% 
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
  }else if(grepl("country", cumulative_type) & grepl("averted_tr", case_type)){
    plt_table <- target_table_country %>% 
      filter(threshold == threshold) %>% 
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
  }else if(grepl("country", cumulative_type) & grepl("averted_ob", case_type)){
    plt_table <- target_table_country %>% 
      filter(threshold == threshold) %>% 
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
  }


  # The third: cumulative histogram by districts 
  if(grepl("district", cumulative_type) & grepl("true", case_type)){
    plt_table <- target_table_mean_district_across_year %>% 
      filter(threshold == threshold) %>%
      filter(general_scenario == "campaign_default" | surveillance_scenario == "district-estimate")
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
      filter(threshold == threshold) %>%
      filter(general_scenario == "campaign_default" | surveillance_scenario == "district-estimate")
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
  }else if(grepl("district", cumulative_type) & grepl("averted_tr", case_type)){
    plt_table <- target_table_mean_district_across_year %>% 
      filter(threshold == threshold)
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
  }else if(grepl("district", cumulative_type) & grepl("averted_ob", case_type)){
    plt_table <- target_table_mean_district_across_year %>% 
      filter(threshold == threshold)
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
  }


  # Return one single plot 
  return(plt)

}



#' @export
#' @name plot_efficacy
#' @title plot_efficacy
#' @description plot the efficacy and return the ggplot object
#' @param cache the cached environment
#' @param compare_type type of the comparison
#' @param threshold threshold to plot  
#' @param cumulative_type either by year (non_cumulative) or cumulated across years
#' @return cached file names 
plot_efficacy <- function(cache, compare_type, threshold, cumulative_type){

  # Restructure of the target table -- no NA's originally because it's all true case
  target_table <- cache$target_table %>%
    mutate(general_scenario = "both", 
          campaign_default_clinical_case = campaign_default_true_case / true_confirm_rate, 
          no_vaccination_clinical_case = no_vaccination_true_case / true_confirm_rate) %>% 
    mutate(averted_true_case = no_vaccination_true_case - campaign_default_true_case, 
          averted_clinical_case = no_vaccination_clinical_case - campaign_default_clinical_case) %>%
    mutate(dose = actual_fvp * 2, 
          efficacy = averted_true_case / dose)
  target_table$efficacy[is.infinite(target_table$efficacy)] <- 0
  target_table$efficacy[is.na(target_table$efficacy)] <- 0
  
  if(grepl("subtract", compare_type)){
    table_no <- target_table[target_table$surveillance_scenario == "no-estimate", ] %>% 
      rename(efficacy_no_estimate = efficacy) %>% 
      select(ISO, admin_level, NAME_1, NAME_2, incid_trend, outbk_trend, threshold, year, run_id, efficacy_no_estimate)
    target_table <- target_table[target_table$surveillance_scenario != "no-estimate", ]
    target_table <- dplyr::left_join(target_table, table_no, by = c("ISO", "admin_level", "NAME_1", "NAME_2", "incid_trend", "outbk_trend", "threshold", "year", "run_id")) %>% 
      mutate(efficacy = efficacy - efficacy_no_estimate) %>% 
      select(- efficacy_no_estimate)
  }


  # For side-by-side comparison and direct comparison 
  if(grepl("side", compare_type) & grepl("non", cumulative_type)){
    target_table_total_country <- target_table %>% 
      group_by(ISO, admin_level, incid_trend, outbk_trend, general_scenario, surveillance_scenario, threshold, year, run_id) %>% 
      summarise(dose = sum(dose), 
                averted_true_case = sum(averted_true_case)) %>%
      mutate(run_id = paste0("total_country", run_id), efficacy = averted_true_case / dose)
    target_table_total_country$efficacy[is.infinite(target_table_total_country$efficacy)] <- 0
    target_table_total_country$efficacy[is.na(target_table_total_country$efficacy)] <- 0
    
  }else if(grepl("side", compare_type) & grepl("year", cumulative_type)){
    target_table_total_country <- target_table %>% 
      group_by(ISO, admin_level, incid_trend, outbk_trend, general_scenario, surveillance_scenario, threshold, run_id) %>% 
      summarise(dose = sum(dose), 
                averted_true_case = sum(averted_true_case)) %>%
      mutate(run_id = paste0("total_country", run_id), efficacy = averted_true_case / dose)
    target_table_total_country$efficacy[is.infinite(target_table_total_country$efficacy)] <- 0
    target_table_total_country$efficacy[is.na(target_table_total_country$efficacy)] <- 0
    
  }else if(grepl("side", compare_type) & grepl("district", cumulative_type)){
    target_table_total_country <- target_table %>% 
      group_by(ISO, admin_level, NAME_1, NAME_2, incid_trend, outbk_trend, general_scenario, surveillance_scenario, threshold, run_id) %>% 
      summarise(dose = sum(dose), 
                averted_true_case = sum(averted_true_case)) %>%
      mutate(efficacy = averted_true_case / dose)
    target_table_total_country$efficacy[is.infinite(target_table_total_country$efficacy)] <- 0
    target_table_total_country$efficacy[is.na(target_table_total_country$efficacy)] <- 0
    
  }else if(grepl("subtract", compare_type) & grepl("non", cumulative_type)){
    target_table_total_country <- target_table %>% 
      group_by(ISO, admin_level, incid_trend, outbk_trend, general_scenario, surveillance_scenario, threshold, year, run_id) %>% 
      summarise(efficacy = weighted.mean(efficacy, dose, na.rm = T)) 
    target_table_total_country$efficacy[is.infinite(target_table_total_country$efficacy)] <- 0
    target_table_total_country$efficacy[is.na(target_table_total_country$efficacy)] <- 0
  
  }else if(grepl("subtract", compare_type) & grepl("year", cumulative_type)){
    target_table_total_country <- target_table %>% 
      group_by(ISO, admin_level, incid_trend, outbk_trend, general_scenario, surveillance_scenario, threshold, run_id) %>% 
      summarise(efficacy = weighted.mean(efficacy, dose, na.rm = T)) 
    target_table_total_country$efficacy[is.infinite(target_table_total_country$efficacy)] <- 0
    target_table_total_country$efficacy[is.na(target_table_total_country$efficacy)] <- 0
  
  }else if(grepl("subtract", compare_type) & grepl("district", cumulative_type)){
    target_table_total_country <- target_table %>% 
      group_by(ISO, admin_level, NAME_1, NAME_2, incid_trend, outbk_trend, general_scenario, surveillance_scenario, threshold, run_id) %>% 
      summarise(efficacy = weighted.mean(efficacy, dose, na.rm = T)) 
    target_table_total_country$efficacy[is.infinite(target_table_total_country$efficacy)] <- 0
    target_table_total_country$efficacy[is.na(target_table_total_country$efficacy)] <- 0
    
  }


  # Plot
  if(grepl("side", compare_type) & grepl("non", cumulative_type)){
    plt_table <- target_table_total_country %>% 
      filter(threshold == threshold) %>% 
      mutate(year = as.factor(year))
    plt <- plt_table %>% 
      ggplot(aes(x=year, y=efficacy, fill=year)) +
      geom_boxplot() + 
      facet_grid( ISO + admin_level ~ surveillance_scenario, scales = "free", space = "free") + 
      theme_minimal() + 
      coord_flip()
    # pdf("/home/kaiyuezou/VIMC_Model/surveillance_project/gavi_vimc_cholera/diagnostics/surveillance_project/test13.pdf")
    # plt
    # dev.off()
  }else if(grepl("side", compare_type) & grepl("year", cumulative_type)){
    plt_table <- target_table_total_country %>% 
      filter(threshold == threshold)
    plt <- plt_table %>% 
      ggplot(aes(x=ISO, y=efficacy, fill=ISO)) +
      geom_boxplot() + 
      facet_grid( ISO + admin_level ~ surveillance_scenario, scales = "free", space = "free") + 
      theme_minimal()
    # pdf("/home/kaiyuezou/VIMC_Model/surveillance_project/gavi_vimc_cholera/diagnostics/surveillance_project/test14.pdf")
    # plt
    # dev.off()
  }else if(grepl("side", compare_type) & grepl("district", cumulative_type)){
    plt_table <- target_table_total_country %>% 
      filter(threshold == threshold)
    plt <- plt_table %>% 
      ggplot(aes(x=NAME_1, y=efficacy, fill=ISO)) +
      geom_boxplot() + 
      facet_grid( ISO + admin_level ~ surveillance_scenario, scales = "free", space = "free") + 
      theme_minimal() + 
      coord_flip()
    # pdf("/home/kaiyuezou/VIMC_Model/surveillance_project/gavi_vimc_cholera/diagnostics/surveillance_project/test17.pdf", height = 20, width = 15)
    # plt
    # dev.off()
  }else if(grepl("subtract", compare_type) & grepl("non", cumulative_type)){
    plt_table <- target_table_total_country %>% 
      filter(threshold == threshold) %>% 
      mutate(year = as.factor(year))
    plt <- plt_table %>% 
      ggplot(aes(x=year, y=efficacy, fill=year)) +
      geom_boxplot() + 
      facet_grid( ISO + admin_level ~ surveillance_scenario, scales = "free", space = "free") + 
      theme_minimal() + 
      coord_flip()
    # pdf("/home/kaiyuezou/VIMC_Model/surveillance_project/gavi_vimc_cholera/diagnostics/surveillance_project/test15.pdf")
    # plt
    # dev.off()
  }else if(grepl("subtract", compare_type) & grepl("year", cumulative_type)){
    plt_table <- target_table_total_country %>% 
      filter(threshold == threshold)
    plt <- plt_table %>% 
      ggplot(aes(x=ISO, y=efficacy, fill=ISO)) +
      geom_boxplot() + 
      facet_grid( ISO + admin_level ~ surveillance_scenario, scales = "free", space = "free") + 
      theme_minimal()
    # pdf("/home/kaiyuezou/VIMC_Model/surveillance_project/gavi_vimc_cholera/diagnostics/surveillance_project/test16.pdf")
    # plt
    # dev.off()
  }else if(grepl("subtract", compare_type) & grepl("district", cumulative_type)){
    plt_table <- target_table_total_country %>% 
      filter(threshold == threshold)
    plt <- plt_table %>% 
      ggplot(aes(x=NAME_1, y=efficacy, fill=ISO)) +
      geom_boxplot() + 
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
