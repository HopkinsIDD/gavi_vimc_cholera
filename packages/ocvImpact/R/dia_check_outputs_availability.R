#' @name check_outputs_availability
#' @title check_outputs_availability
#' @description function to check if all the output files are there 
#' @param rawoutpath path of output_raw in which all the model output files are stored
#' @param configpath path of configuration files
#' @param countries vector of all the countries included the diagnostic report
#' @param scenario takes value of either "campaign-default" or "no-vaccination"
#' @param threshold threshold chosen for the diagnostic report, takes one of 1/1000, 1/5000, or 1/10,000
#' @param surveillance_scenario scale of surveillance, takes "no-estimate", "global-estimate", or "district-estimate"
#' @return TRUE or FALSE indicating whether the output files for the specified countries and scenarios exist in output folder, and also return warning messages to check output of which country is missing
#' @export
#' @include
check_outputs_availability <- function(rawoutpath = "output_raw/202110gavi-3",
                                        configpath = "configs/202110gavi-3",
                                        countries,
                                        scenario,
                                        surveillance_scenario,
                                        threshold,
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
    
    if(!file.exists(filename_admin1) | !file.exists(filename_admin2)){
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

# test check_outputs_availability()
# check_outputs_availability(countries = c("NGA"), scenario = "campaign-default", surveillance_scenario = "district-estimate", threshold = 0.001)
