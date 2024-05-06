
#' @title prepare_config
#' @name prepare_config
#' @description Automate generation of model configuration files from a dataframe of parameters. Generates one config per row in the p dataframe.
#'
#' @param p dataframe with config parameters, one config per row
#' @param configpath path to config files
#' @return Message returning name of yaml config file 
#' @export
prepare_config <- function(p, configpath){
  if(!p$targeting == "threshold_unconstrained" & !p$runname == "202310gavi-4"){
    config_name <- paste0(configpath, "/", paste(p$country, p$scenario, p$nsamples, sep = "_"), ".yml")
    sink(file = config_name)

    cat(paste0(
      "runname: '", p$runname, "'\n",
      "country: '", p$country, "'\n",
      "scenario: '", p$scenario, "'\n",
      "clean: ", p$clean, "\n",
      "incid:\n",
      "  num_samples: ", p$nsamples, "\n",
      "  redraw: ", p$redrawIncid, "\n",
      "  use_country_incid_trend: ", p$use_country_incid_trend, "\n",
      "vacc:\n",
      "  targeting_strategy: ", p$targeting, "\n",
      "  num_skip_years: ", p$nskipyear, "\n",
      "setting:\n",
      "  incidence_rate_trend: ", p$incidence_rate_trend, "\n",
      "  outbreak_multiplier: ", p$outbreak_multiplier, "\n", 
      "  random_seed: ", p$random_seed, "\n"
    ))

    sink()

  } else if(!p$targeting == "threshold_unconstrained" & p$runname == "202310gavi-4"){
    ## for the DRC Case Study
     if (p$use_montagu_coverage == FALSE){
      
       config_name <- paste0(configpath, "/", paste(p$country, p$scenario, p$nsamples, p$ndoses, p$use_montagu_coverage, p$targeting, sep = "_"), ".yml")
       sink(file = config_name)
      
       cat(paste0(
        "runname: '", p$runname, "'\n",
        "country: '", p$country, "'\n",
        "scenario: '", p$scenario, "'\n",
        "clean: ", p$clean, "\n",
        "campaign_cov: ", p$campaign_cov, "\n",
        "incid:\n",
        "  num_samples: ", p$nsamples, "\n",
        "  redraw: ", p$redrawIncid, "\n",
        "  use_country_incid_trend: ", p$use_country_incid_trend, "\n",
        "vacc:\n",
        "  targeting_strategy: ", p$targeting, "\n",
        "  num_skip_years: ", p$nskipyear, "\n",
        "  ndoses: ", p$ndoses, "\n",
        "custom:\n",
        "  use_montagu_coverage: ", p$use_montagu_coverage, "\n",
        "  output_years: ", p$output_years, "\n",
        "  targeting_filename: ", p$custom_targeting_filename, "\n",
        "  coverage_filename: ", p$custom_coverage_filename, "\n",
        "  use_custom_shapefile: ", p$use_custom_shapefile, "\n",
        "  shapefile_filename: ", p$custom_shapefile_filename, "\n",
        "  country_shapefile_filename: ", p$custom_country_shapefile_filename, "\n",
        "setting:\n",
        "  incidence_rate_trend: ", p$incidence_rate_trend, "\n",
        "  outbreak_multiplier: ", p$outbreak_multiplier, "\n", 
        "  random_seed: ", p$random_seed, "\n"
       ))
      
       sink()  
     } else {
      config_name <- paste0(configpath, "/", paste(p$country, p$scenario, p$nsamples, p$ndoses, sep = "_"), ".yml")
      sink(file = config_name)
     
      cat(paste0(
       "runname: '", p$runname, "'\n",
       "country: '", p$country, "'\n",
       "scenario: '", p$scenario, "'\n",
       "clean: ", p$clean, "\n",
       "incid:\n",
       "  num_samples: ", p$nsamples, "\n",
       "  redraw: ", p$redrawIncid, "\n",
       "  use_country_incid_trend: ", p$use_country_incid_trend, "\n",
       "vacc:\n",
       "  targeting_strategy: ", p$targeting, "\n",
       "  num_skip_years: ", p$nskipyear, "\n",
       "  ndoses: ", p$ndoses, "\n",
       "setting:\n",
       "  incidence_rate_trend: ", p$incidence_rate_trend, "\n",
       "  outbreak_multiplier: ", p$outbreak_multiplier, "\n", 
       "  random_seed: ", p$random_seed, "\n"
     ))
     
     sink()
     
     }
    
  }else if(p$targeting == "threshold_unconstrained"){
    config_name <- paste0(configpath, "/", paste(p$country, p$scenario, p$surveillance_scenario, p$nsamples, sep = "_"), ".yml") #for now 
    sink(file = config_name)

    cat(paste0(
      "runname: '", p$runname, "'\n",
      "country: '", p$country, "'\n",
      "scenario: '", p$scenario, "'\n",
      "clean: ", p$clean, "\n",
      "incid:\n",
      "  num_samples: ", p$nsamples, "\n",
      "  redraw: ", p$redrawIncid, "\n",
      "  use_country_incid_trend: ", p$use_country_incid_trend, "\n",
      "vacc:\n",
      "  targeting_strategy: '", p$targeting, "'\n",
      "  vac_incid_threshold: ", p$vac_incid_threshold, "\n",
      "  vac_unconstrained: ", p$vac_unconstrained, "\n",
      "  vac_admin_level: '", p$vac_admin_level, "'\n",
      "  vac_coverage: ", p$vac_coverage, "\n",
      "  num_skip_years: ", p$nskipyear, "\n",
      "  ndoses: ", p$ndoses, "\n",
      "  vac_interval: ", p$vac_interval, "\n",
      "  sim_start_year: ", p$sim_start_year, "\n",
      "  vac_start_year: ", p$vac_start_year, "\n",
      "  vac_end_year: ", p$vac_end_year, "\n",
      "  sim_end_year: ", p$sim_end_year, "\n",
      "  use_mean_ir: ", p$use_mean_ir, "\n",
      "  mean_ir_span: ", p$mean_ir_span, "\n",
      "surveillance_scenario:\n",
      "  surveillance_scenario: '", p$surveillance_scenario, "'\n",
      "  testing_sensitivity: ", p$testing_sensitivity, "\n",
      "setting:\n",
      "  incidence_rate_trend: ", p$incidence_rate_trend, "\n",
      "  outbreak_multiplier: ", p$outbreak_multiplier, "\n", 
      "  random_seed: ", p$random_seed, "\n", 
      "optimize:\n",
      "  save_intermediate_raster: ", p$save_intermediate_raster, "\n", 
      "  save_final_output_raster: ", p$save_final_output_raster, "\n", 
      "  ir_pre_screening: ", p$ir_pre_screening, "\n"
    ))

    sink()
  }

  return(message(paste(config_name, "written")))

}


