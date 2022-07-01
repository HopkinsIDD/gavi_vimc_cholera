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
                                threshold, 
                                ...){   
  
  # make a table for the output
  df_setting <- matrix(NA, nrow = length(countries), ncol = 13) %>% as.data.frame()
  names(df_setting) <- c("country", "scenario", "surveillance_scenario", "threshold", "num_samples", "vac_coverage", "num_skip_years", "sim_start_year", "vac_start_year", "vac_end_year", "sim_end_year", "incidence_trend", "outbreak_multiplier")
  df_setting$country <- countries

  # check if the desiring config files are available in config folder
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
  return(df_setting)
}

# test check_model_setting()
# check_model_setting(countries = c("NGA"), scenario = "campaign-default", surveillance_scenario = "district-estimate", threshold = 0.001)


