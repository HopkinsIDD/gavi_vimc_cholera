#' @name run_surveillance_scenario
#' @title run_surveillance_scenario
#' @description Run full surveillance project scenario
#' @param datapath path to input data 
#' @param modelpath path to montagu files
#' @param country country code
#' @param scenario Unique string that identifies the coverage scenario name
#' @param rawoutpath path to raw model output files
#' @param nsamples number of stochastic samples to use
#' @param ve_direct vaccine effect function (default: generate_pct_protect_function())
#' @param indirect_mult indirect effects function (default: generate_indirect_incidence_mult())
#' @param secular_trend_mult function that returns a multiplier for secular incidence trends over time (default: generate_flatline_multiplier())
#' @param clean logical that indicates whether existing targeting and model outputs (sus, pop, vacc) should be deleted (default = TRUE)
#' @param redraw logical that indicates whether existing incidence raster samples should be redrawn (default = FALSE)
#' @param config the whole config file that will be read in and used 
#' @param ... Optional parameters to pass to [`assign_vaccine_targets()`]. See [`assign_vaccine_targets()`] for defaults.
#' @return 
#' @export
#' @include utils.R create_static_modelInputs.R create_sus_modelInputs.R create_expectedCases.R
run_surveillance_scenario <- function(
  datapath, 
  modelpath, 
  country, 
  scenario,
  rawoutpath,
  nsamples, 
  ve_direct = generate_pct_protect_function(),
  indirect_mult = generate_indirect_incidence_mult(),
  secular_trend_mult = function(a,b,c,d){return(a*b*c*d)},
  clean = TRUE,
  redraw = FALSE,
  config, 
  ...){

  ##### Get all the parameters ready 
  save_intermediate_raster <- as.logical(config$optimize$save_intermediate_raster)
  save_final_output_raster <- as.logical(config$optimize$save_final_output_raster)
  targeting_strategy <- config$vacc$targeting_strategy
  vac_incid_threshold <- as.numeric(config$vacc$vac_incid_threshold)
  vac_unconstrained <- as.logical(config$vacc$vac_unconstrained)
  vac_admin_level <- tolower(config$vacc$vac_admin_level)
  vac_coverage <- as.numeric(config$vacc$vac_coverage)
  surveillance_scenario <- config$surveillance_scenario$surveillance_scenario
  vac_interval <- as.numeric(config$vacc$vac_interval)
  sim_start_year <- as.numeric(config$vacc$sim_start_year)
  vac_start_year <- as.numeric(config$vacc$vac_start_year)
  vac_end_year <- as.numeric(config$vacc$vac_end_year)
  sim_end_year <- as.numeric(config$vacc$sim_end_year)
  num_skip_years <- as.numeric(config$vacc$num_skip_years)

  incidence_rate_trend <- as.logical(config$setting$incidence_rate_trend)
  use_country_incid_trend <- as.logical(config$incid$use_country_incid_trend)
  outbreak_multiplier <- as.logical(config$setting$outbreak_multiplier)  
  random_seed <- as.numeric(config$setting$random_seed)



  ##### Initialize the campaign table and load the incidence rate raster 
  if(vac_admin_level != "both"){stop("The only supported vaccination campagin administration level is 'both', please check the config. ")}
  shp0 <- load_shp0_by_country(datapath, country)
  shp1 <- load_shp1_by_country(datapath, country)
  shp2 <- load_shp2_by_country(datapath, country)
  rc_list <- load_baseline_incidence(datapath, modelpath, country, campaign_cov = vac_coverage, baseline_year = sim_start_year, first_vacc_year = vac_start_year, 
                                     incidence_rate_trend, use_country_incid_trend, shp0 = shp0, shp1 = shp1, shp2 = shp2, 
                                     random_seed = random_seed, nsamples = nsamples, redraw = redraw)



  ##### Loop through all the years to simulate 
  for(model_year in sim_start_year:sim_end_year){
    
    #### Update the list of districts to target at the beginning of a year -- this is where stochasticity introduced due to random draw of confirmation rates
    message(paste("Starting simulation of year", model_year, "in country", country))
    for(layer_idx in 1:nsamples){
      rc_list[[layer_idx]] <- update_targets_list(datapath = datapath, modelpath = modelpath, country = country, scenario = scenario, 
                                                  rc_list = rc_list[[layer_idx]], model_year = model_year, campaign_cov = vac_coverage, 
                                                  threshold = vac_incid_threshold, vac_unconstrained = vac_unconstrained, 
                                                  surveillance_scenario = surveillance_scenario, 
                                                  vac_interval = vac_interval, 
                                                  vac_start_year = vac_start_year, vac_end_year = vac_end_year, 
                                                  num_skip_years = num_skip_years)
    }


    #### Update the vaccinated proportion raster 
    ### Read in pop raster for this year
    pop <- ocvImpact::create_model_pop_raster(datapath, modelpath, country, model_year)
    
    ### Calculate/update pop and vacc raster input
    if(!exists("input_list") | save_intermediate_raster){
      input_list_exp <- parse(text = paste0( "input_list <- list(", toString(rep("NULL", nsamples)), ")" ))
      eval(input_list_exp)
    }

    for(layer_idx in 1:nsamples){
      input_list[[layer_idx]] <- update_vac_raster( datapath, modelpath, country, scenario, rawoutpath,
                                                    rc_list = rc_list[[layer_idx]], model_year = model_year, pop = pop,
                                                    input_list = input_list[[layer_idx]])
    } 
    
    if(save_intermediate_raster){  
      save_vac_raster(datapath, modelpath, country, nsamples, model_year, input_list)
      input_list <- NULL
    }


    #### Calculate/update suspectible population raster
    #### Now only save the newest year 
    if(!save_intermediate_raster){
      sus_list_exp <- parse(text = paste0( "sus_list <- list(", toString(rep("NULL", nsamples)), ")" ))
      eval(sus_list_exp)

      for(layer_idx in 1:nsamples){
        sus_list[[layer_idx]] <- update_sus_rasterStack(datapath, modelpath, country, scenario, rawoutpath,
                                                        pop = pop,
                                                        ve_direct = ve_direct,
                                                        baseline_year = sim_start_year,
                                                        model_year = model_year, 
                                                        input_list = input_list[[layer_idx]], # generated from update_input_rasterStack() function
                                                        sus_list = sus_list[[layer_idx]] # rasterStack of proportion susceptible generated in last year
                                                        )
      }
    }else{
      sus_list <- update_sus_rasterStack_optimized( datapath, modelpath, country, scenario, rawoutpath,
                                                    pop = pop,
                                                    ve_direct = ve_direct,
                                                    baseline_year = sim_start_year,
                                                    model_year = model_year
                                                  )
      save_sus_raster(datapath, modelpath, country, nsamples, model_year, sus_list)
      sus_list <- NULL
    }


    #### Get the expected cases for the year
    # indirect_mult = generate_indirect_incidence_mult() #temp, will be deleted **************************************
    # secular_trend_mult = function(a,b,c,d){return(a*b*c*d)} #temp, will be deleted **************************************
    # if(!exists("ec_list")){ec_list <- NULL}
    ec_list <- surveillance_create_expectedCases( datapath, 
                                                  modelpath, 
                                                  country, 
                                                  scenario, 
                                                  rawoutpath, 
                                                  vacc_alloc = NULL, 
                                                  indirect_mult, 
                                                  secular_trend_mult, 
                                                  nsamples, 
                                                  is_cf = TRUE, 
                                                  clean = clean, 
                                                  redraw, 
                                                  sus_list, 
                                                  pop, 
                                                  model_year,
                                                  config)


    #### Get ready for the next year -- this needs to be optimized 
    if(model_year < sim_end_year){
      for(layer_idx in 1:nsamples){
        ec_list_tmp <- list("ec_rasterStack_admin1" = ec_list$ec_rasterStack_admin1[[layer_idx]],
                            "ec_rasterStack_admin2" = ec_list$ec_rasterStack_admin2[[layer_idx]])
        rc_list[[layer_idx]] <- surveillance_add_rc_new_row(rc_list[[layer_idx]], ec_list = ec_list_tmp, pop, model_year, sim_start_year, sim_end_year, shp1, shp2)
        rm(ec_list_tmp)
      }
    }
    rm(ec_list)
    
  }



  ##### Write output files
  if(save_final_output_raster){
  # vaccination proportion and susceptible proportion rasterStack
    for(model_year in sim_start_year:sim_end_year){
      save_vac_raster(datapath, modelpath, country, nsamples, model_year, input_list, 
                      rawoutpath, clean, scenario, incidence_rate_trend, outbreak_multiplier, vac_incid_threshold, surveillance_scenario)
      save_sus_raster(datapath, modelpath, country, nsamples, model_year, sus_list, 
                      rawoutpath, clean, scenario, incidence_rate_trend, outbreak_multiplier, vac_incid_threshold, surveillance_scenario)
    }
  }

  # # proportion susceptible rasterStack
  # dir.create(paste0(rawoutpath, "/", scenario, "/"), showWarnings = FALSE)
  # sus1_out_fn <- paste0(rawoutpath, "/", scenario, "/", paste("incid", incidence_rate_trend, "outbk", outbreak_multiplier, 
  #                       vac_incid_threshold, surveillance_scenario, country, sep = "_"), "_sus_admin1.tif")
  # sus2_out_fn <- paste0(rawoutpath, "/", scenario, "/", paste("incid", incidence_rate_trend, "outbk", outbreak_multiplier, 
  #                       vac_incid_threshold, surveillance_scenario, country, sep = "_"), "_sus_admin2.tif")
  # message(paste("Writing proportion susceptible rasterStack for", country))
  # if( !file.exists(sus1_out_fn) | (file.exists(sus1_out_fn)&clean) ){
  #   raster::writeRaster(sus_list[["sus_rasterStack_admin1"]], filename = sus1_out_fn, overwrite = TRUE)}
  # if( !file.exists(sus2_out_fn) | (file.exists(sus2_out_fn)&clean) ){
  #   raster::writeRaster(sus_list[["sus_rasterStack_admin2"]], filename = sus2_out_fn, overwrite = TRUE)}
  
  # # expected cases rasterStack
  # ec1_out_fn <- paste0( rawoutpath, "/", scenario, "/", paste("incid", incidence_rate_trend, "outbk", outbreak_multiplier, 
  #                       vac_incid_threshold, surveillance_scenario, country, sep = "_"), "_ec_admin1.tif")
  # ec2_out_fn <- paste0( rawoutpath, "/", scenario, "/", paste("incid", incidence_rate_trend, "outbk", outbreak_multiplier, 
  #                       vac_incid_threshold, surveillance_scenario, country, sep = "_"), "_ec_admin2.tif")
  # message(paste("Writing expected cases rasterStack for", country))
  # if( !file.exists(ec1_out_fn) | (file.exists(ec1_out_fn)&clean) ){
  #   raster::writeRaster(ec_list[["ec_rasterStack_admin1"]], filename = ec1_out_fn, overwrite = TRUE)}
  # if( !file.exists(ec2_out_fn) | (file.exists(ec2_out_fn)&clean) ){
  #   raster::writeRaster(ec_list[["ec_rasterStack_admin2"]], filename = ec2_out_fn, overwrite = TRUE)}
  
  # rc_list with incidence and targeted area of each modeled year -- only save the first one for now
  rc1_out_fn <- paste0(rawoutpath, "/", scenario, "/", paste("incid", incidence_rate_trend, "outbk", outbreak_multiplier, 
                        vac_incid_threshold, surveillance_scenario, country, sep = "_"), "_rc_admin1.csv")
  rc2_out_fn <- paste0(rawoutpath, "/", scenario, "/", paste("incid", incidence_rate_trend, "outbk", outbreak_multiplier, 
                        vac_incid_threshold, surveillance_scenario, country, sep = "_"), "_rc_admin2.csv")
  rc1 <- rc_list[[1]]$rc1 %>% sf::st_drop_geometry()
  rc2 <- rc_list[[1]]$rc2 %>% sf::st_drop_geometry()
  message(paste("Writing targeting tables of all years for", country))
  if( !file.exists(rc1_out_fn) | (file.exists(rc1_out_fn)&clean) ){
    readr::write_csv(rc1, rc1_out_fn)}
  if( !file.exists(rc2_out_fn) | (file.exists(rc2_out_fn)&clean) ){
    readr::write_csv(rc2, rc2_out_fn)}
  
  
  
  message(paste("End of simulating vaccination campaigns in", country, "from", sim_start_year, "to", sim_end_year))
  rm(rc_list, input_list, sus_list, ec_list, shp1, shp2)
  gc()
  
  return(NULL)

}


