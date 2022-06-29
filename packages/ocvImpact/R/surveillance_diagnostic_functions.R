###################
# Functions
###################

get_filenames <- function(cache, surveillance_project_directory, pre_country, pre_vac_incid_thresholds, 
                          no_vaccination_surveillance_scenario){
  # Get the set_all_parameters.R to have the context first
  source(file.path(surveillance_project_directory, "scripts/set_all_parameters.R"))

  # Get the correct version of country list and threshold list
  country_list <- ifelse(pre_country %in% countries, pre_country, countries)  
  rm(pre_country, countries)

  threshold_list <- ifelse(pre_vac_incid_thresholds %in% vac_incid_thresholds, pre_vac_incid_thresholds, vac_incid_thresholds)
  rm(pre_vac_incid_thresholds, vac_incid_thresholds)

  admin_list <- ifelse(vac_admin_level == "both", c("admin1", "admin2"), vac_admin_level)
  year_list <- sim_start_year:sim_end_year

  # Generate and cache the file names 
  ec_filenames_camp <- apply(expand.grid(paste0("output_raw/", runname, "/campaign-default/", "incid"), 
                            incidence_rate_trend, "outbk", outbreak_multiplier, threshold_list, surveillance_scenarios, 
                            country_list[1], "ec", admin_list, paste0(year_list, ".tif")), 1, paste, collapse = "_")
  ec_filenames_novac <- apply(expand.grid(paste0("output_raw/", runname, "/no-vaccination/", "incid"), 
                            incidence_rate_trend, "outbk", outbreak_multiplier, threshold_list, no_vaccination_surveillance_scenario, 
                            country_list[1], "ec", admin_list, paste0(year_list, ".tif")), 1, paste, collapse = "_")
  ecraster_fns <- c(ec_filenames_camp, ec_filenames_novac)

  ectable_fns_camp <- apply(expand.grid(paste0("output_raw/", runname, "/campaign-default/", country_list[1]), 
                            "ec", paste0(admin_list, ".csv")), 1, paste, collapse = "_")
  ectable_fns_novac <- apply(expand.grid(paste0("output_raw/", runname, "/no-vaccination/", country_list[1]), 
                            "ec", paste0(admin_list, ".csv")), 1, paste, collapse = "_")                          
  ectable_fns <- c(ectable_fns_camp, ectable_fns_novac)

  ## this needs to be changed! or not
  targettable_fns_camp <- apply(expand.grid(paste0("output_raw/", runname, "/campaign-default/", "incid"), 
                            incidence_rate_trend, "outbk", outbreak_multiplier, threshold_list, surveillance_scenarios, 
                            country_list[1], "rc", admin_list, paste0(year_list, ".csv")), 1, paste, collapse = "_")
  targettable_fns_novac <- apply(expand.grid(paste0("output_raw/", runname, "/no-vaccination/", "incid"), 
                            incidence_rate_trend, "outbk", outbreak_multiplier, threshold_list, no_vaccination_surveillance_scenario, 
                            country_list[1], "rc", admin_list, paste0(year_list, ".csv")), 1, paste, collapse = "_")
  targettable_fns <- c(targettable_fns_camp, targettable_fns_novac)

  truecomfirmrate_fns <- paste0("intermediate_raster/", country_list, "_trueconfirmrate_admin2.tif")                          

}
