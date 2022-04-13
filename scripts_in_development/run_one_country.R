# This function runs simulation for one country
run_one_country <- function(datapath, country, campaign_cov = 0.8, modelpath, rawoutpath,
                            baseline_year = 2018, end_year, nsamples = 30, redraw = FALSE){
  
  # load country shp at both levels
  shp1 <- load_shp1_by_country(datapath, country)
  shp2 <- load_shp2_by_country(datapath, country)
  
  # load baseline incidence within the template table
  rc_list <- load_baseline_incidence(datapath, country, campaign_cov, baseline_year,
                                     shp1 = shp1, shp2 = shp2)
  
  # prepare for entering loop
  lambda <- create_incid_raster(modelpath, datapath, country, nsamples = 30, redraw = FALSE) # load baseline incidence rasterStack
  model_year <- baseline_year
  
  ## start the loop: one loop for each campaign year
  for(model_year in baseline_year:end_year){
    
    message(paste("Starting simulation of campaign of year", model_year))
    # select the targets
    rc_list <- update_targets_list(rc_list = rc_list, baseline_year = baseline_year, campaign_cov = campaign_cov)
    
    # read in pop raster for this year
    latest_campaign_year <- max(rc_list[[1]]$target_year)
    pop <- create_model_pop_raster(datapath, modelpath, country, latest_campaign_year) # might need to pre-run function to get this outside this function, because in the functions creating population rasters (in order to update incidence), this will be reused
    
    
    # calculate/update pop and vacc raster input
    if(latest_campaign_year == baseline_year){
      
      input_list <- create_first_year_input(datapath, modelpath, country,
                                            rc_list = rc_list, baseline_year = baseline_year, pop = pop)
    } else{
      
      input_list <- update_input_rasterStack(datapath,modelpath, country, scenario, rawoutpath,
                                             rc_list = rc_list, pop = pop,
                                             baseline_year = baseline_year, this_year = latest_campaign_year,
                                             input_list = input_list) # input an empty list for the first year, for the following year, input is that list from last year)
    }
    
    
    # calculate/update suspectible raster 
    ve_direct <- generate_pct_protect_function()
    if(latest_campaign_year == baseline_year){
      
      sus_list <- create_first_year_sus(datapath, modelpath, country, scenario, rawoutpath,
                                        rc_list = rc_list, ve_direct = ve_direct,
                                        baseline_year = baseline_year, 
                                        pop = pop,
                                        input_list = input_list, # generated from update_input_rasterStack() function
                                        clean)
      
    } else{
      
      sus_list <- update_sus_rasterStack(datapath, modelpath, country, scenario, rawoutpath,
                                         rc_list = rc_list,
                                         pop = pop,
                                         ve_direct = ve_direct,
                                         baseline_year = baseline_year,
                                         this_year = latest_campaign_year, 
                                         input_list = input_list, # generated from update_input_rasterStack() function
                                         sus_list = sus_list  # rasterStack of proportion susceptible generated in last year
                                         )
    }
    
    
    # get expected case raster
    indirect_mult <- generate_indirect_incidence_mult()
    secular_trend_mult <- generate_flatline_multiplier()
    
    if(latest_campaign_year == baseline_year){
      ec_list <- create_first_year_ec(datapath, modelpath, country, scenario, rawoutpath,
                                      rc_list = rc_list, 
                                      baseline_year = baseline_year,
                                      sus_list = sus_list, # susceptible rasterStack generated in the last function
                                      input_list = input_list,      # population and proportion vaccinated rasterStack generated in this loop
                                      indirect_mult =indirect_mult,
                                      secular_trend_mult = secular_trend_mult,
                                      lambda = lambda # loaded baseline incidence of this country, using: lambda <- create_incid_raster(modelpath, datapath, country, nsamples, redraw)
                                      # nsamples, 
                                      # is_cf,
                                      # redraw
                                      )
    } else{
      ec_list <- update_ec(datapath, modelpath, country, scenario, rawoutpath,
                           rc_list = rc_list, 
                           baseline_year = baseline_year,
                           sus_list = sus_list, # susceptible rasterStack generated in the last function
                           input_list = input_list,      # population and proportion vaccinated rasterStack generated in this loop
                           indirect_mult =indirect_mult,
                           secular_trend_mult = secular_trend_mult,
                           lambda = lambda, # loaded baseline incidence of this country, using: lambda <- create_incid_raster(modelpath, datapath, country, nsamples, redraw)
                           # nsamples, 
                           # is_cf,
                           # redraw,
                           ec_list = ec_list)
    }
    
    
    # append empty rows to prepare for next campaign
    rc_list <- new_campaign_preparation(rc_list = rc_list, shp1 = shp1, shp2 = shp2)
    
    # fill in new incidence for the new campaign year 
    rc_list <- update_incidence(rc_list, ec_list)
    #write.csv(rc_list[[1]] %>% st_drop_geometry(), "~/vimc/draft/output/rc_list_admin1_loop.csv")
    #write.csv(rc_list[[2]] %>% st_drop_geometry(), "~/vimc/draft/output/rc_list_admin2_loop.csv")
    
    message(paste("Campaign of year", latest_campaign_year, "completed, incidence of year", latest_campaign_year+1, "updated and ready for new round of targeting."))
  
    model_year <- model_year + 1
  } # for loop end (campaign of all years end)
  
  
  ## Write output files
  
  # proportion susceptible rasterStack
  sus1_out_fn <- paste0(rawoutpath, "/", scenario, "/", country, "_sus_admin1.tif")
  sus2_out_fn <- paste0(rawoutpath, "/", scenario, "/", country, "_sus_admin2.tif")
  message(paste("Writing proportion susceptible rasterStack for", country))
  raster::writeRaster(sus_list[["sus_rasterStack_admin1"]], filename = sus1_out_fn)
  raster::writeRaster(sus_list[["sus_rasterStack_admin2"]], filename = sus2_out_fn)
  
  
  # expected cases rasterStack
  ec1_out_fn <- paste0(rawoutpath, "/", scenario, "/", country, "_ec_admin1.tif")
  ec2_out_fn <- paste0(rawoutpath, "/", scenario, "/", country, "_ec_admin2.tif")
  message(paste("Writing expected cases rasterStack for", country))
  raster::writeRaster(ec_list[["ec_rasterStack_admin1"]], filename = ec1_out_fn)
  raster::writeRaster(ec_list[["ec_rasterStack_admin2"]], filename = ec2_out_fn)
  
  
  # rc_list with incidence and targeted area of each modeled year
  rc1_out_fn <- paste0(rawoutpath, "/", scenario, "/", country, "_rc_admin1.csv")
  rc2_out_fn <- paste0(rawoutpath, "/", scenario, "/", country, "_rc_admin2.csv")
  rc1 <- rc_list[["rc1"]] %>% filter(target_year <= end_year) %>% st_drop_geometry()
  rc2 <- rc_list[["rc2"]] %>% filter(target_year <= end_year) %>% st_drop_geometry()
  message(paste("Writing targeting tables of all years for", country))
  write.csv(rc1, rc1_out_fn)
  write.csv(rc2, rc2_out_fn)
  
  
  message(paste("End of simulating vaccination campaigns in", country, "from", baseline_year, "to", end_year))
  rm(rc_list, input_list, sus_list, ec_list, shp1, shp2, lambda)
  gc()
  
  return(NULL)
  
}
