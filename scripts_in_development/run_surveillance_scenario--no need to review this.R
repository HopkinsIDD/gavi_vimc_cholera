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



  ##### Initialize the campaign table and load the incidence rate raster 
  if(vac_admin_level != "both"){stop("The only supported vaccination campagin administration level is 'both', please check the config. ")}
  shp1 <- load_shp1_by_country(datapath, country)
  shp2 <- load_shp2_by_country(datapath, country)
  rc_list <- load_baseline_incidence(datapath, country, campaign_cov = vac_coverage, baseline_year = sim_start_year, first_vacc_year = vac_start_year, 
                                     shp1 = shp1, shp2 = shp2)
  
  # lambda <- ocvImpact::create_incid_raster(modelpath, datapath, country, nsamples, redraw)



  ##### Loop through all the years to simulate 
  for(model_year in sim_start_year:sim_end_year){
    
    #### Update the list of districts to target at the beginning of a year -- this is where stochasticity introduced due to random draw of confirmation rates
    message(paste("Starting simulation of campaign of year", model_year, "in country", country))
    rc_list <- update_targets_list( datapath = datapath, modelpath = modelpath, country = country, scenario = scenario, 
                                    rc_list = rc_list, model_year = model_year, baseline_year = vac_start_year, campaign_cov = vac_coverage, 
                                    threshold = vac_incid_threshold, vac_unconstrained = vac_unconstrained, 
                                    surveillance_scenario = surveillance_scenario, 
                                    vac_interval = vac_interval, 
                                    vac_start_year = vac_start_year, vac_end_year = vac_end_year, 
                                    num_skip_years = num_skip_years)


    #### Update the vaccinated proportion raster 
    ### Read in pop raster for this year
    # latest_campaign_year <- max(rc_list[[1]]$target_year)
    pop <- ocvImpact::create_model_pop_raster(datapath, modelpath, country, model_year)
    
    ### Calculate/update pop and vacc raster input
    if(model_year == sim_start_year){
      
      input_list <- create_first_year_vac_raster( datapath, modelpath, country,
                                                  rc_list = rc_list, model_year = model_year, pop = pop)
    }else{
      
      input_list <- update_vac_raster(datapath,modelpath, country, scenario, rawoutpath,
                                      rc_list = rc_list, pop = pop,
                                      model_year,
                                      input_list = input_list) # input an empty list for the first year, for the following year, input is that list from last year)
    }
    

    #### Calculate/update suspectible population raster 
    ve_direct <- generate_pct_protect_function() #for temp use
      
    if(!exists("sus_list")){sus_list <- NULL}
    sus_list <- update_sus_rasterStack( datapath, modelpath, country, scenario, rawoutpath,
                                        rc_list = rc_list,
                                        pop = pop,
                                        ve_direct = ve_direct,
                                        baseline_year = sim_start_year,
                                        model_year = model_year, 
                                        input_list = input_list, # generated from update_input_rasterStack() function
                                        sus_list = sus_list # rasterStack of proportion susceptible generated in last year
                                      )
    

    #### Get the expected cases for the year
    indirect_mult = generate_indirect_incidence_mult() #temp, will be deleted 
    secular_trend_mult = function(a,b,c,d){return(a*b*c*d)} #temp, will be deleted 
    if(!exists("ec_list")){ec_list <- NULL}
    expCases <- surveillance_create_expectedCases(datapath, 
                                                  modelpath, 
                                                  country, 
                                                  scenario, 
                                                  rawoutpath, 
                                                  vacc_alloc = NULL, 
                                                  indirect_mult, 
                                                  secular_trend_mult, 
                                                  nsamples, 
                                                  is_cf = TRUE, 
                                                  redraw, 
                                                  sus_list, 
                                                  pop, 
                                                  model_year, 
                                                  ec_list = ec_list, 
                                                  config)
    
    

    

  }


  vacc_alloc <- allocate_vaccine(datapath, modelpath, country, scenario, ...) #the changes start from here 

  ## write proportion vaccinated to file and export total population raster stack
  dummy <- create_static_modelInputs(datapath, modelpath, country, scenario, rawoutpath, vacc_alloc, clean)

  ## write susceptible population proportion raster 
  dummy2 <- create_sus_modelInputs(datapath, modelpath, country, scenario, rawoutpath, vacc_alloc, ve_direct, clean)

  incidence_rate_trend <- as.logical(config$setting$incidence_rate_trend)
  outbreak_multiplier <- as.logical(config$setting$outbreak_multiplier)
  setting <- paste0('incid_trend_', incidence_rate_trend, '_outb_layer_',  outbreak_multiplier)

  dir.create(paste0(rawoutpath, "/", scenario, "/", setting), showWarnings = FALSE)
  ec_out_fn <- paste0(rawoutpath, "/", scenario, "/", setting, "/", country, "_ec.csv")
  if(clean | !file.exists(ec_out_fn)){ ## rerun

    if (is.null(vacc_alloc)){
      message("Calculate expected cases: with vaccination")
      expCases <- create_expectedCases(datapath, modelpath, country, scenario, rawoutpath, vacc_alloc, indirect_mult, secular_trend_mult, nsamples, is_cf = TRUE, redraw)
    } else{
      message("Calculate expected cases: no vaccination")
      expCases <- create_expectedCases(datapath, modelpath, country, scenario, rawoutpath, vacc_alloc, indirect_mult, secular_trend_mult, nsamples, is_cf = FALSE, redraw)
    } 

      ## Write to file 
      message(paste("Write expected cases:", country, scenario, "\n", ec_out_fn))
      readr::write_csv(expCases, ec_out_fn)

  } else{ ## read existing
    message(paste("Reading expected cases:", country, scenario, "\n", ec_out_fn))
    expCases <- readr::read_csv(ec_out_fn)
  }
  



  return(expCases)

}



