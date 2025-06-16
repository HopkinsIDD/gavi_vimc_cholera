#' @name surveillance_create_expectedCases
#' @title surveillance_create_expectedCases
#' @description Create proportion of population vaccinated and total population rasterStacks. Write vaccination raster to file and export population raster
#' @param datapath path to input data 
#' @param modelpath path to montagu files
#' @param country country code
#' @param scenario Unique string that identifies the coverage scenario name
#' @param rawoutpath path to raw model output files
#' @param indirect_mult function for calculating indirect effects
#' @param secular_trend_mult function for calculating secular trends in mean annual incidence
#' @param nsamples number of stochastic samples to use
#' @param clean
#' @param redraw logical indicate whether to redraw incidence raster samples; pass to [`create_incid_raster()`]
#' @param sus_list
#' @param pop
#' @param model_year
#' @param config
#' @param rc_targeted
#' @return
#' @export
#' @include utils.R utils_montagu.R create_incid_raster.R load_shapefile_by_country.R 
surveillance_create_expectedCases <- function(
  datapath,
  modelpath,
  country,
  scenario,
  rawoutpath,
  indirect_mult,
  secular_trend_mult,
  nsamples,
  clean, 
  redraw, 
  sus_list, 
  pop, 
  model_year,
  config, 
  rc_targeted
  ){

  message("Calculating expected cases now. ")
  message(Sys.time())
  ############# Use the configs to determine the setting #############
  incidence_rate_trend <- as.logical(config$setting$incidence_rate_trend)
  use_country_incid_trend <- as.logical(config$incid$use_country_incid_trend)
  outbreak_multiplier <- as.logical(config$setting$outbreak_multiplier)
  sim_start_year <- as.numeric(config$vacc$sim_start_year)
  sim_end_year <- as.numeric(config$vacc$sim_end_year)

  vac_incid_threshold <- as.numeric(config$vacc$vac_incid_threshold)
  surveillance_scenario <- config$surveillance_scenario$surveillance_scenario
  save_intermediate_raster <- as.logical(config$optimize$save_intermediate_raster)
  save_final_output_raster <- as.logical(config$optimize$save_final_output_raster)
  
  
  ### Get the rasters and shapefile ready 
  lambda <- create_incid_raster(modelpath, datapath, country, nsamples, redraw = redraw, random_seed = config$setting$random_seed)
  pop_rasterLayer <- pop
  rm(pop)
  shp0 <- load_shapefile_by_country(datapath, country, simple = TRUE)
  
  if(!save_intermediate_raster){
    sus_rasterLayer1 <- sus_list[[1]]$sus_rasterStack_admin1
    if(!is.null(sus_rasterLayer1)){
      if(nsamples != 1){
        for(layer_idx in 2:nsamples){
          sus_rasterLayer1 <- c(sus_rasterLayer1, sus_list[[layer_idx]]$sus_rasterStack_admin1)
        }        
      }
    }
    sus_rasterLayer2 <- sus_list[[1]]$sus_rasterStack_admin2
    if(!is.null(sus_rasterLayer2)){
      if(nsamples != 1){
        for(layer_idx in 2:nsamples){
          sus_rasterLayer2 <- c(sus_rasterLayer2, sus_list[[layer_idx]]$sus_rasterStack_admin2)
        } 
      }
    }
    rm(sus_list)

  }else{
    if(scenario == "campaign-default"){
      sus_admin1_fn <- paste0("intermediate_raster/", country, "_sus_admin1_", model_year, ".tif")
      sus_admin2_fn <- paste0("intermediate_raster/", country, "_sus_admin2_", model_year, ".tif")
    
      if("rc1" %in% rc_targeted){sus_rasterLayer1 <- terra::rast(sus_admin1_fn)}else{sus_rasterLayer1 <- NULL}
      if("rc2" %in% rc_targeted){sus_rasterLayer2 <- terra::rast(sus_admin2_fn)}else{sus_rasterLayer2 <- NULL}
    }else{
      sus_rasterLayer1 <- sus_list$sus_rasterStack_admin1 
      sus_rasterLayer2 <- sus_list$sus_rasterStack_admin2
    }

  }
  
  # #fixing the MRT issue
  # if(country == 'MRT'){
  #   lambda <- raster::setExtent(lambda, raster::extent(shp0), keepres=FALSE, snap=FALSE)
  #   pop_rasterLayer <- raster::setExtent(pop_rasterLayer, raster::extent(shp0), keepres=FALSE, snap=FALSE)
  # }
  # #fixing the BGD issue
  # if(country == 'BGD'){
  #   lambda <- raster::resample(lambda, pop_rasterLayer, method = "ngb")
  # }
  

  ### The incidence rate trend multiplier 
  if(incidence_rate_trend){
    incid_trend_function <- ocvImpact::generate_flatline_multiplier(
                                        trendtype = 'incidence rate', 
                                        datapath = datapath, 
                                        modelpath = modelpath, 
                                        country = country, 
                                        use_country_incid_trend = use_country_incid_trend)
  }else{ 
    incid_trend_function <- function(year){return(1)}
  }
  
 
  ### The outbreak multiplier and the model year 
  ## make sure that the directory where all outbreak data is saved exists 
  dir.create(paste0(datapath, '/outbreak'), showWarnings = FALSE)
  oy <- model_year

  ## get the multiplier function 
  if(outbreak_multiplier){
    # first check if the outbreak raster file already exists, generate one if not
    random_seed <- as.numeric(config$setting$random_seed)
    setting_num <- random_seed #for now, yes
    outbreak_out_fn <- paste0(datapath, '/outbreak/', country, '_', (oy %% 10), '.tif')

    if(!file.exists(outbreak_out_fn)){
      outbreak_trend_function <- ocvImpact::outbreak_incidence_rate_multiplier(
                                            datapath = datapath,
                                            modelpath = modelpath,
                                            country = country,
                                            lambda = lambda, 
                                            population_raster = pop_rasterStack, 
                                            output_years = c(oy), 
                                            use_country_incid_trend = use_country_incid_trend)
    }else{
      outbreak_multiplier_raster <- c(outbreak_out_fn) #please note that the raster here is RasterStack, not RasterBrick
      outbreak_trend_function <- function(yr_index){return(outbreak_multiplier_raster)}
    }
    
  }else{
    outbreak_trend_function <- function(yr_index){return(1)}
  }


  ### Get the overall multiplier ready 
  incid_trend_multiplier <- incid_trend_function(year = model_year)
  outbreak_ic_multiplier <- outbreak_trend_function(yr_index = 1) ### the returned list only has one element to use anyways
  confirmation_multiplier <- surveillance_true_confirmation_rate(country, admin_level = 2)
  utilization_multiplier <- 1 ### 11/25 added
  
  overall_multiplier <- secular_trend_mult(incid_trend_multiplier, outbreak_ic_multiplier, confirmation_multiplier, utilization_multiplier)
    

  ### Multiply all the layers 

  ## admin1 first 
  # make new indirect effects template
  if("rc1" %in% rc_targeted){
    indirect_rasterLayer <- sus_rasterLayer1
    terra::values(indirect_rasterLayer) <- indirect_mult(1-as.numeric(terra::values(sus_rasterLayer1)))

    ec_rasterStack1 <- tryCatch(
      if(!is.numeric(overall_multiplier) & inherits(overall_multiplier, "SpatRaster")){
        lambda * sus_rasterLayer1 * indirect_rasterLayer * pop_rasterLayer * overall_multiplier # SpatRaster could be multiplied directly
      }else{
        lambda * sus_rasterLayer1 * indirect_rasterLayer * pop_rasterLayer * overall_multiplier
      },

      error = function(e){
        print(paste0('The year when it goes wrong is ', oy))
      }
    )
    ec_rasterStack1 <- c(ec_rasterStack1)
  }else{ec_rasterStack1 <- NULL}

  ## admin2 then 
  # make new indirect effects template
  if("rc2" %in% rc_targeted){
    indirect_rasterLayer <- sus_rasterLayer2
    terra::values(indirect_rasterLayer) <- indirect_mult(1-as.numeric(terra::values(sus_rasterLayer2)))

    ec_rasterStack2 <- tryCatch(
      if(!is.numeric(overall_multiplier) & inherits(overall_multiplier, "SpatRaster")){
        lambda * sus_rasterLayer2 * indirect_rasterLayer * pop_rasterLayer * overall_multiplier # SpatRaster could be multiplied directly
      }else{
        lambda * sus_rasterLayer2 * indirect_rasterLayer * pop_rasterLayer * overall_multiplier
      },

      error = function(e){
        print(paste0('The year when it goes wrong is ', oy))
      }
    )
    ec_rasterStack2 <- c(ec_rasterStack2)
  }else{ec_rasterStack2 <- NULL}


  ### Save the raster -- true cases ********************************************
  if(save_final_output_raster){
    ec1_out_fn <- paste0(rawoutpath, "/", scenario, "/", paste("incid", incidence_rate_trend, "outbk", outbreak_multiplier, 
                        vac_incid_threshold, surveillance_scenario, country, sep = "_"), "_ec_admin1_", model_year, ".tif")
    ec2_out_fn <- paste0(rawoutpath, "/", scenario, "/", paste("incid", incidence_rate_trend, "outbk", outbreak_multiplier, 
                        vac_incid_threshold, surveillance_scenario, country, sep = "_"), "_ec_admin2_", model_year, ".tif")
    
    message(paste("Writing expected true cases rasterStack for", country))
    dir.create(paste0(rawoutpath, "/", scenario, "/"), showWarnings = FALSE)
    if( !file.exists(ec1_out_fn) | (file.exists(ec1_out_fn)&clean) ){
      if(!is.null(ec_rasterStack1)){terra::writeRaster(ec_rasterStack1, filename = ec1_out_fn, overwrite = TRUE)}}
    if( !file.exists(ec2_out_fn) | (file.exists(ec2_out_fn)&clean) ){
      if(!is.null(ec_rasterStack2)){terra::writeRaster(ec_rasterStack2, filename = ec2_out_fn, overwrite = TRUE)}}
  
  }

  
  ### optimize memory usage
  rm(outbreak_ic_multiplier)
  rm(outbreak_trend_function)
  rm(sus_rasterLayer1)
  rm(sus_rasterLayer2)
  rm(lambda)
  rm(indirect_rasterLayer)
  rm(overall_multiplier)
  gc()


  # ### Save the tables -- true cases here ********************************************
  # ec_yr1 <- exactextractr::exact_extract(ec_rasterStack1, shp0, fun = "sum", stack_apply = TRUE)
  # ec_yr2 <- exactextractr::exact_extract(ec_rasterStack2, shp0, fun = "sum", stack_apply = TRUE)
  # ec_vec1 <- as.numeric(ec_yr1)
  # ec_vec2 <- as.numeric(ec_yr2)
  # pop_total <- as.numeric(exactextractr::exact_extract(pop_rasterLayer, shp0, fun = "sum", stack_apply = TRUE))
  # mean_incid1 <- ec_vec1 / pop_total
  # mean_incid2 <- ec_vec2 / pop_total

  # mean_incid_vec1 <- as.numeric(mean_incid1)
  # mean_incid_vec2 <- as.numeric(mean_incid2)
  # rc1 <- tibble::tibble(country = country, incidence_rate_trend = incidence_rate_trend, outbreak_multiplier = outbreak_multiplier, 
  #                       vac_incid_threshold = vac_incid_threshold, surveillance_scenario = surveillance_scenario, 
  #                       year = oy, run_id = seq_along(ec_vec1), ec = ec_vec1, incid_rate = mean_incid_vec1)
  # rc2 <- tibble::tibble(country = country, incidence_rate_trend = incidence_rate_trend, outbreak_multiplier = outbreak_multiplier, 
  #                       vac_incid_threshold = vac_incid_threshold, surveillance_scenario = surveillance_scenario, 
  #                       year = oy, run_id = seq_along(ec_vec2), ec = ec_vec2, incid_rate = mean_incid_vec2)

  # ## External dataset 
  # dir.create(paste0(rawoutpath, "/", scenario, "/"), showWarnings = FALSE)
  # ec_out_fn1 <- paste0(rawoutpath, "/", scenario, "/", country, "_ec_admin1.csv")
  # ec_out_fn2 <- paste0(rawoutpath, "/", scenario, "/", country, "_ec_admin2.csv")
  
  # if(file.exists(ec_out_fn1)){
  #   ec_out1 <- readr::read_csv(ec_out_fn1)
  #   if(nrow(ec_out1[ec_out1$country == country & ec_out1$incidence_rate_trend == incidence_rate_trend 
  #                   & ec_out1$outbreak_multiplier == outbreak_multiplier 
  #                   & ec_out1$vac_incid_threshold == vac_incid_threshold 
  #                   & ec_out1$surveillance_scenario == surveillance_scenario 
  #                   & ec_out1$year == oy, ])>0 & clean){
  #     ec_out1 <- ec_out1[!(ec_out1$country == country & ec_out1$incidence_rate_trend == incidence_rate_trend 
  #                         & ec_out1$outbreak_multiplier == outbreak_multiplier 
  #                         & ec_out1$vac_incid_threshold == vac_incid_threshold 
  #                         & ec_out1$surveillance_scenario == surveillance_scenario 
  #                         & ec_out1$year == oy), ] 
  #   }
  #   ec_out1 <- rbind(ec_out1, rc1)
    
  # }else{
  #   ec_out1 <- rc1
  # }

  # if(file.exists(ec_out_fn2)){
  #   ec_out2 <- readr::read_csv(ec_out_fn2)
  #   if(nrow(ec_out2[ec_out2$country == country & ec_out2$incidence_rate_trend == incidence_rate_trend 
  #                   & ec_out2$outbreak_multiplier == outbreak_multiplier 
  #                   & ec_out2$vac_incid_threshold == vac_incid_threshold 
  #                   & ec_out2$surveillance_scenario == surveillance_scenario 
  #                   & ec_out2$year == oy, ])>0 & clean){
  #     ec_out2 <- ec_out2[!(ec_out2$country == country & ec_out2$incidence_rate_trend == incidence_rate_trend 
  #                         & ec_out2$outbreak_multiplier == outbreak_multiplier 
  #                         & ec_out2$vac_incid_threshold == vac_incid_threshold 
  #                         & ec_out2$surveillance_scenario == surveillance_scenario 
  #                         & ec_out2$year == oy), ]
  #   }
  #   ec_out2 <- rbind(ec_out2, rc2)

  # }else{
  #   ec_out2 <- rc2
  # }

  # readr::write_csv(ec_out1, ec_out_fn1)
  # readr::write_csv(ec_out2, ec_out_fn2)
  

  ### Return the list -- true cases ********************************************
  message("Finished calculating expected cases and generating the new ec_list that contains suspected cases. ")
  message(Sys.time())
  ec_list <- list("ec_rasterStack_admin1" = ec_rasterStack1,
                  "ec_rasterStack_admin2" = ec_rasterStack2)
  rm(ec_rasterStack1)
  rm(ec_rasterStack2)

  return(ec_list)

}


