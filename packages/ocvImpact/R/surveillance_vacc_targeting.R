#' @name load_shp0_by_country
#' @title load_shp0_by_country
#' @description load_shp0_by_country
#' @param datapath path to input data 
#' @param country country code 
#' @return 
#' @export
#' @include
load_shp0_by_country <- function(datapath, country){
  
  tryCatch(
    {
      country_pattern <- paste(country, "0", sep = "_")
      shp <- GADMTools::gadm_sf_loadCountries(c(country), level = 0, basefile = file.path(datapath, "shapefiles/"))$sf
      
      message(paste0("Loading ", datapath, "/shapefiles/", country_pattern, "_sf.rds"))
    },
    error = function(e) {
      print(paste("Unable to get shapefile for", country, ".", e))
    }
  )
  
  return(shp)
}

#' @name load_shp1_by_country
#' @title load_shp1_by_country
#' @description load_shp1_by_country
#' @param datapath path to input data 
#' @param country country code 
#' @param simple return country level shapefile is TRUE  
#' @return 
#' @export
#' @include
load_shp1_by_country <- function(datapath, country, simple = FALSE){
  
  tryCatch(
    {
      if (simple){
        country_pattern <- paste(country, "0", sep = "_")
        shp <- GADMTools::gadm_sf_loadCountries(c(country), level = 0, basefile = file.path(datapath, "shapefiles/"))$sf
      } else{
        country_pattern <- paste(country, "1", sep = "_")
        shp <- GADMTools::gadm_sf_loadCountries(c(country), level = 1, basefile = file.path(datapath, "shapefiles/"))$sf
      }
      message(paste0("Loading ", datapath, "/shapefiles/", country_pattern, "_sf.rds"))
    },
    error = function(e) {
      print(paste("Unable to get shapefile for", country, ".", e))
    }
  )
  
  return(shp)
}

#' @name load_shp2_by_country
#' @title load_shp2_by_country
#' @description load_shp2_by_country
#' @param datapath path to input data 
#' @param country country code 
#' @param simple return country level shapefile is TRUE  
#' @return 
#' @export
#' @include
load_shp2_by_country <- function(datapath, country, simple = FALSE){
  
  tryCatch(
    {
      if (simple){
        country_pattern <- paste(country, "0", sep = "_")
        shp <- GADMTools::gadm_sf_loadCountries(c(country), level = 0, basefile = file.path(datapath, "shapefiles/"))$sf
      } else{
        country_pattern <- paste(country, "2", sep = "_")
        shp <- GADMTools::gadm_sf_loadCountries(c(country), level = 2, basefile = file.path(datapath, "shapefiles/"))$sf
      }
      message(paste0("Loading ", datapath, "/shapefiles/", country_pattern, "_sf.rds"))
    },
    error = function(e) {
      print(paste("Unable to get shapefile for", country, ".", e))
    }
  )
  
  return(shp)
}

#' @name load_baseline_incidence
#' @title load_baseline_incidence
#' @description from lambda to incidence rate look-up table to rc_list 
#' @param datapath
#' @param modelpath
#' @param country  
#' @param campaign_cov
#' @param baseline_year
#' @param first_vacc_year
#' @param incidence_rate_trend
#' @param use_country_incid_trend
#' @param shp0
#' @param shp1
#' @param shp2
#' @param random_seed
#' @param nsamples
#' @param redraw
#' @return 
#' @export
#' @include
load_baseline_incidence <- function(datapath, 
                                    modelpath, 
                                    country, 
                                    campaign_cov, 
                                    baseline_year, 
                                    first_vacc_year, 
                                    incidence_rate_trend, 
                                    use_country_incid_trend, 
                                    shp0, 
                                    shp1, 
                                    shp2, 
                                    random_seed, 
                                    nsamples, 
                                    redraw){
  
  ## load the incidence rate raster and initialize the rc_list 
  country_baseline <- create_incid_raster(modelpath, datapath, country, nsamples, redraw, random_seed)
  rc_list <- NULL

  ## load population data of the first year 
  pop_baseline <- ocvImpact::create_model_pop_raster(datapath, modelpath, country, baseline_year)
  
  ## mask the incidence raster to population raster first 
  country_baseline <- raster::mask(country_baseline, pop_baseline)
  pop_baseline <- raster::mask(pop_baseline, country_baseline)[[1]]

  ## if the incidence rate trend should be implemented 
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
  country_baseline <- country_baseline*incid_trend_function(baseline_year)
  
  ## get the pop ready
  pop1 <- get_admin_population(pop_baseline, shp1)
  pop2 <- get_admin_population(pop_baseline, shp2)
  # total population
  total_pop <- sum(pop1)
  # total population test
  if(sum(pop1) != sum(pop2)){warning(paste(" *** population summarized across admin1 and admin2 districts differ by", abs(sum(pop1) - sum(pop2))))}
  
  ## get the trur confirm rate layer ready 
  dir.create(paste0("intermediate_raster/"), showWarnings = FALSE)
  confirm_rate_fn <- paste0("intermediate_raster/", country, "_trueconfirmrate_admin2.tif")
  
  if(file.exists(confirm_rate_fn)){
    message("The true confirm rate rasters have already been in place, they will not be replaced. ")
    confirm_rate_raster <- raster::stack(confirm_rate_fn)
    if(raster::nlayers(confirm_rate_raster) != nsamples){
      stop("The number of layers of the true confirm rate raster in the intermediate folder is wrong, run killed, please delete the wrong raster. ")
    }
  }else{
    message("Generating the new true confirm rate rasters and saving them after. ")
  }
    
  omicron_dataset <- readr::read_csv(paste0(datapath, "/confirmation_rate/parameters.csv"))
  raster1_template <- raster::calc(country_baseline[[1]], fun = function(x){ifelse(!is.na(x), 1, NA)})
  # prepare to generate the true confirm rate and set the random seed here to make sure they are consistent across different runs/scenarios 
  set.seed(random_seed) #******************************************************************************** IMPORTANT ********************************************************************************
  

  for(layer_idx in 1:nsamples){
    
    ### get the right distribution 
    # confirm_rate_value <- MCMCglmm::rtnorm(n = nrow(shp2), mean = omicron_dataset$mean, sd = omicron_dataset$sd, lower = 0, upper = 1) #truncated normal distribution
    # confirm_rate_value <- runif(n = nrow(shp2), min = 0, max = 1) #uniform distribution 
    # confirm_rate_value <- rnorm(n = nrow(shp2), mean = omicron_dataset$mean, sd = omicron_dataset$sd) #discarded function
    # updated 8/30/22 using the predictive distribution
    confirm_rate_value <- rbeta(n = nrow(shp2), shape1 = omicron_dataset$shape1, shape2 = omicron_dataset$shape2) 


    ### get clinical incid2 first to get the true confirm rate raster and clinical case raster 
    country_baseline_single_layer <- country_baseline[[layer_idx]]
    incid2 <- pop_weighted_admin_mean_incid(datapath, modelpath, incidence_rate_raster = country_baseline_single_layer, pop_raster = pop_baseline, country, admin_shp = shp2)

    confirm_rate_df <- dplyr::mutate(shp2, 
                                    true_confirm_rate = confirm_rate_value, 
                                    # true_incidence_rate = incid2 * true_confirm_rate, 
                                    # pop_model = pop2,
                                    # confirmed_incidence_rate = incid2, 
                                    observed_case = incid2 * pop2) 
    
    ## rasterize
    if(!file.exists(confirm_rate_fn)){
      confirm_rate_raster <- fasterize::fasterize(
          confirm_rate_df,
          raster1_template,
          field = "true_confirm_rate",
          fun = "last",
          background = 1
      )
      confirm_rate_raster <- raster::mask(confirm_rate_raster, raster1_template, updatevalue = NA)
      confirm_rate_raster <- raster::mask(confirm_rate_raster, pop_baseline, updatevalue = NA)
      confirm_rate_raster_one_layer <- confirm_rate_raster
    }else{
      confirm_rate_raster_one_layer <- confirm_rate_raster[[layer_idx]]
    }
    
    observed_case_raster <- fasterize::fasterize(
        confirm_rate_df,
        raster1_template,
        field = "observed_case",
        fun = "last",
        background = 1
    )
    observed_case_raster <- raster::mask(observed_case_raster, raster1_template, updatevalue = NA)
    confirm_rate_value_for_admin1 <- pop_weighted_admin_mean_incid(datapath, modelpath, confirm_rate_raster_one_layer, observed_case_raster, country, admin_shp = shp1)
    

    ### Get the true incidence rates
    country_baseline_single_layer_true <- country_baseline_single_layer * confirm_rate_raster_one_layer
    incid1_true <- pop_weighted_admin_mean_incid(datapath, modelpath, incidence_rate_raster = country_baseline_single_layer_true, pop_raster = pop_baseline, country, admin_shp = shp1)
    incid2_true <- pop_weighted_admin_mean_incid(datapath, modelpath, incidence_rate_raster = country_baseline_single_layer_true, pop_raster = pop_baseline, country, admin_shp = shp2)
    ## clean up
    rm(confirm_rate_raster_one_layer, observed_case_raster, country_baseline_single_layer, confirm_rate_raster_one_layer, country_baseline_single_layer_true)


    ### Make the df
    # get incidence and pop table (admin 1 level)
    rc1 <- dplyr::mutate(shp1, 
                        true_incidence_rate = incid1_true,
                        confirmation_lens = NA,
                        confirmation_rate = NA, 
                        confirmed_incidence_rate = NA,  
                        mean_confirmed_incidence_rate = NA,  
                        pop_model = pop1,
                        pop_prop = pop1/total_pop,
                        year = baseline_year, 
                        latest_target_year = NA,
                        is_target = NA,  # whether that place was target in this year
                        actual_prop_vaccinated = NA,
                        actual_fvp = NA, 
                        true_confirm_rate = confirm_rate_value_for_admin1) %>% # fully vaccinated population
      dplyr::mutate(  true_incidence_rate = ifelse(is.na(true_incidence_rate), 0, true_incidence_rate), 
                      true_case = NA, 
                      mean_true_incidence_rate = true_incidence_rate ) %>% 
      dplyr::select(ISO, NAME_0, NAME_1, true_confirm_rate, true_incidence_rate, mean_true_incidence_rate, 
                    confirmation_lens, confirmation_rate, confirmed_incidence_rate, mean_confirmed_incidence_rate, 
                    pop_model, pop_prop, year, latest_target_year, is_target,
                    actual_prop_vaccinated, actual_fvp, true_case)
    
    # get incidence and pop table (admin 2 level)
    rc2 <- dplyr::mutate(shp2, 
                        true_incidence_rate = incid2_true,
                        confirmation_lens = NA,
                        confirmation_rate = NA, 
                        confirmed_incidence_rate = NA,  
                        mean_confirmed_incidence_rate = NA,  
                        pop_model = pop2,
                        pop_prop = pop2/total_pop,
                        year = baseline_year, 
                        latest_target_year = NA,
                        is_target = NA,  # whether that place was target in this year
                        actual_prop_vaccinated = NA,
                        actual_fvp = NA, 
                        true_confirm_rate = confirm_rate_value) %>% # fully vaccinated population
      dplyr::mutate(  true_incidence_rate = ifelse(is.na(true_incidence_rate), 0, true_incidence_rate), 
                      true_case = NA, 
                      mean_true_incidence_rate = true_incidence_rate ) %>% 
      dplyr::select(ISO, NAME_0, NAME_1, NAME_2, true_confirm_rate, true_incidence_rate, mean_true_incidence_rate, 
                    confirmation_lens, confirmation_rate, confirmed_incidence_rate, mean_confirmed_incidence_rate, 
                    pop_model, pop_prop, year, latest_target_year, is_target,
                    actual_prop_vaccinated, actual_fvp, true_case)
    
    # Save the confirm_rate_raster raster into the intermediate folder
    if(!file.exists(confirm_rate_fn)){
      if(layer_idx == 1){
        tmp <- confirm_rate_raster
      }else if(layer_idx < nsamples){
        tmp <- raster::stack(tmp, confirm_rate_raster)
      }else{
        tmp <- raster::stack(tmp, confirm_rate_raster)
        if(!file.exists(confirm_rate_fn)){ raster::writeRaster(tmp, filename = confirm_rate_fn, overwrite = FALSE) }
        rm(tmp, confirm_rate_raster)
      }
    }

    
    # check on population proportion 
    if (sum(rc1$pop_prop)!=1 | sum(rc2$pop_prop)!=1){
      warning(paste(" *** The population proportion calculation is incorrect for", country, ", they are", sum(rc1$pop_prop), "and", sum(rc2$pop_prop)))
    }
    
    # save the tables 
    rc_list_single_layer <- list("rc1" = rc1, "rc2" = rc2)
    if(is.null(rc_list)){
      rc_list <- list(rc_list_single_layer)
    }else{
      rc_list[[length(rc_list)+1]] <- rc_list_single_layer
    }

  }
  

  ## Clean up and return 
  rm(rc_list_single_layer)
  rm(pop1, pop2)
  rm(rc1, rc2)
  gc()
  
  return(rc_list)
}

#' @name update_targets_list
#' @title update_targets_list
#' @description update targets in the rc_list
#' @param datapath path to input data 
#' @param modelpath 
#' @param country
#' @param scenario 
#' @param rc_list 
#' @param model_year
#' @param campaign_cov
#' @param threshold
#' @param vac_unconstrained
#' @param surveillance_scenario
#' @param vac_interval country level skipped years 
#' @param vac_start_year
#' @param vac_end_year 
#' @param num_skip_years district level skipped years 
#' @param rc_targeted
#' @param use_mean_ir
#' @param mean_ir_span
#' @param testing_sensitivity true sensitivity of testing 
#' @return 
#' @export
#' @include
update_targets_list <- function(datapath, modelpath, country, scenario, 
                                rc_list, model_year, 
                                campaign_cov, 
                                threshold, vac_unconstrained, 
                                surveillance_scenario, 
                                vac_interval, #this is country level 
                                vac_start_year, vac_end_year, 
                                num_skip_years, #this is district level 
                                rc_targeted, 
                                use_mean_ir,  #whether use the mean true ir averaged across years 
                                mean_ir_span,  #how many years to average upon 
                                testing_sensitivity
                                ){
  # Decide if just target one admin level 
  if(length(rc_targeted) == 2){
    ### Assertain that this table has new empty rows to fill in and it's the model year 
    if(max(rc_list$rc1$year) != model_year | max(rc_list$rc2$year) != model_year){
      stop("Cannot update the rc list because the year with empty row does not agree with the model year. ")
    }
    
    ### What is the vacc strategy? 
    if(!vac_unconstrained){
      stop("The constrained surveillance scenario is not developed yet. ") #future: maybe some years to skip the vacc(constrained in this sense)? 
    }

    ### What if this is no vac scenario or the vac_interval is larger than 1 year or the vacc is skipped for the year or the vacc has ended or not yet started?
    if(
        ( (!all(is.na(rc_list$rc1$latest_target_year)) & !all(is.na(rc_list$rc2$latest_target_year)))
            & (model_year - max(rc_list$rc1$latest_target_year, na.rm = TRUE) < vac_interval 
                & model_year - max(rc_list$rc2$latest_target_year, na.rm = TRUE) < vac_interval) )      | 
        (scenario == "no-vaccination")                                                                  |
        (model_year > vac_end_year)                                                                     |
        (model_year < vac_start_year) ){
      
      message(paste("The vaccination campaign is skipped for year", model_year, "for the whole country", country))
      
      if( (max(rc_list$rc1$latest_target_year, na.rm = TRUE) != max(rc_list$rc2$latest_target_year, na.rm = TRUE)) & (model_year <= vac_end_year) &
          (   (model_year - max(rc_list$rc1$latest_target_year, na.rm = TRUE) >= vac_interval) 
            | (model_year - max(rc_list$rc2$latest_target_year, na.rm = TRUE) >= vac_interval)   ) ){
        message("admin1 and admin2 do not agree on the latest target year in the rc list table, please check. ")
        cat(paste("The country level vaccination elapse year for the admin1 scenario is", (model_year - max(rc_list$rc1$latest_target_year, na.rm = TRUE))))
        cat(paste("The country level vaccination elapse year for the admin2 scenario is", (model_year - max(rc_list$rc2$latest_target_year, na.rm = TRUE))))
        stop("Error: inconsistency between admin1 and admin2 regarding whether or not to carry out country-level campaign. ")
      }
      rc_list <- update_novacc_year(datapath, rc_list, model_year, surveillance_scenario = surveillance_scenario, rc_targeted, use_mean_ir, mean_ir_span, testing_sensitivity)
      
    }else{
      rc_list <- update_vacc_year(datapath, modelpath, country, rc_list, model_year, 
                                  campaign_cov, threshold, surveillance_scenario, 
                                  vac_start_year, vac_end_year, num_skip_years, rc_targeted, 
                                  use_mean_ir, mean_ir_span, testing_sensitivity)
    }

  # if just one admin level 
  }else{
    if(max(rc_list[[rc_targeted]]$year) != model_year){
      stop("Cannot update the rc list because the year with empty row does not agree with the model year. ")
    }
    if(!vac_unconstrained){
      stop("The constrained surveillance scenario is not developed yet. ") #future: maybe some years to skip the vacc(constrained in this sense)? 
    }

    if( ( (!all(is.na(rc_list[[rc_targeted]]$latest_target_year)))
            & (model_year - max(rc_list[[rc_targeted]]$latest_target_year, na.rm = TRUE) < vac_interval) )         | 
        (scenario == "no-vaccination")                                                                  |
        (model_year > vac_end_year)                                                                     |
        (model_year < vac_start_year) ){
      
      message(paste("The vaccination campaign is skipped for year", model_year, "for the whole country", country))
      rc_list <- update_novacc_year(datapath, rc_list, model_year, surveillance_scenario = surveillance_scenario, rc_targeted, use_mean_ir, mean_ir_span, testing_sensitivity)
      
    }else{
      rc_list <- update_vacc_year(datapath, modelpath, country, rc_list, model_year, 
                                  campaign_cov, threshold, surveillance_scenario, 
                                  vac_start_year, vac_end_year, num_skip_years, rc_targeted, 
                                  use_mean_ir, mean_ir_span, testing_sensitivity)
    }
  }
  
  
  return(rc_list) # updated with latest targets of this year
}

#' @name update_novacc_year
#' @title update_novacc_year
#' @description update the targets in the list if it's not a vaccination campaign year 
#' @param datapath
#' @param rc_list
#' @param model_year
#' @param surveillance_scenario
#' @param rc_targeted
#' @param use_mean_ir
#' @param mean_ir_span
#' @param testing_sensitivity
#' @return 
#' @export
#' @include
update_novacc_year <- function(datapath, rc_list, model_year, surveillance_scenario, rc_targeted, use_mean_ir, mean_ir_span, testing_sensitivity){
  for(rcs in rc_targeted){
    # Update the latest target year, if all NA's from last year or this is the first simulation year, no need to update
    if(model_year > min(rc_list[[rcs]]$year)){ 
      if(!all(is.na(rc_list[[rcs]][rc_list[[rcs]]$year == model_year-1, ]$latest_target_year))){
        rc_list[[rcs]][rc_list[[rcs]]$year == model_year, ]$latest_target_year <- rc_list[[rcs]][rc_list[[rcs]]$year == model_year-1, ]$latest_target_year
      }
    }
    # Update the other variables 
    omicron_dataset <- readr::read_csv(paste0(datapath, "/confirmation_rate/parameters.csv")) #the true confirmation rate 
    rc_list <- get_confirmed_incidence_rate(rc_list, model_year, surveillance_scenario, omicron_dataset, rcs, use_mean_ir, mean_ir_span, testing_sensitivity)

    rc_list[[rcs]][rc_list[[rcs]]$year == model_year, ]$confirmation_lens <- surveillance_scenario
    rc_list[[rcs]][rc_list[[rcs]]$year == model_year, c("is_target", "actual_prop_vaccinated", "actual_fvp")] <- 0 
  }
  return(rc_list)
}

#' @name update_vacc_year
#' @title update_vacc_year
#' @description update the targets in the list if it's a vaccination campaign year 
#' @param datapath
#' @param modelpath
#' @param country
#' @param rc_list
#' @param model_year
#' @param campaign_cov
#' @param threshold
#' @param surveillance_scenario
#' @param vac_start_year
#' @param vac_end_year
#' @param num_skip_years
#' @param rc_targeted
#' @param use_mean_ir
#' @param mean_ir_span
#' @param testing_sensitivity
#' @return 
#' @export
#' @include
update_vacc_year <- function( datapath, modelpath, country, rc_list, model_year, 
                              campaign_cov, threshold, surveillance_scenario, 
                              vac_start_year, vac_end_year, num_skip_years, rc_targeted, 
                              use_mean_ir, mean_ir_span, testing_sensitivity){
  ### Only targeting the population above 1 year old
  if(!"pop_by_age" %in% names(cache)){
    source(paste0("scripts/montagu_handle.R"))
    pop_by_age <- ocvImpact::import_country_agePop(modelpath, country, redownload = FALSE)
    cache$pop_by_age <- pop_by_age
  }else{
    pop_by_age <- cache$pop_by_age
  }
  pop_by_age_1 <- pop_by_age %>% 
    filter(year == model_year) %>% 
    mutate(above_1 = case_when( age_from >= 1 ~ TRUE, 
                                T ~ FALSE)) %>%
    group_by(above_1) %>% 
    summarise(prop = sum(prop_age))
  above_1_proportion <- pop_by_age_1[pop_by_age_1$above_1, ]$prop
  
  ### Put on the confirmation lens 
  ## Read in
  message(paste0("Loading ", datapath, "/confirmation_rate/parameters.csv"))
  omicron_dataset <- readr::read_csv(paste0(datapath, "/confirmation_rate/parameters.csv")) #the true confirmation rate 
  
  ### Use the mean of the confirmed incidence rate or not?
  if(use_mean_ir){targeting_factor <- "mean_confirmed_incidence_rate"}else{targeting_factor <- "confirmed_incidence_rate"}

  ### Loop through the targeted admin levels 
  for (rcs in rc_targeted){
    ## Assign
    rc_list <- get_confirmed_incidence_rate(rc_list, model_year, surveillance_scenario, omicron_dataset, rcs, use_mean_ir, mean_ir_span, testing_sensitivity)

    ### Update the rest of the variables 
    ## Update the target and vaccinated pop
    #if there is no previous year
    if(min(rc_list[[rcs]]$year) == model_year){
      rc_list[[rcs]][rc_list[[rcs]]$year == model_year, ]$is_target <- ifelse(rc_list[[rcs]][rc_list[[rcs]]$year == model_year, ][[targeting_factor]] >= threshold, 1, 0)
    #if there is previous year 
    }else{
      rc_list[[rcs]][rc_list[[rcs]]$year == model_year, ]$is_target <- ifelse(  (!is.na(rc_list[[rcs]][rc_list[[rcs]]$year == model_year-1, ]$latest_target_year)
                                                                                  & model_year - rc_list[[rcs]][rc_list[[rcs]]$year == model_year-1, ]$latest_target_year >= num_skip_years
                                                                                  & rc_list[[rcs]][rc_list[[rcs]]$year == model_year, ][[targeting_factor]] >= threshold)
                                                                              | (is.na(rc_list[[rcs]][rc_list[[rcs]]$year == model_year-1, ]$latest_target_year)
                                                                                  & rc_list[[rcs]][rc_list[[rcs]]$year == model_year, ][[targeting_factor]] >= threshold), 1, 0)
    }
    rc_list[[rcs]][rc_list[[rcs]]$year == model_year, ]$actual_prop_vaccinated <- rc_list[[rcs]][rc_list[[rcs]]$year == model_year, ]$is_target * campaign_cov * above_1_proportion
    rc_list[[rcs]][rc_list[[rcs]]$year == model_year, ]$actual_fvp <- rc_list[[rcs]][rc_list[[rcs]]$year == model_year, ]$actual_prop_vaccinated * rc_list[[rcs]][rc_list[[rcs]]$year == model_year, ]$pop_model
    
    ## Update the latest target year
    if(min(rc_list[[rcs]]$year) == model_year){
      rc_list[[rcs]][rc_list[[rcs]]$year == model_year, ]$latest_target_year <- ifelse(rc_list[[rcs]][rc_list[[rcs]]$year == model_year, ]$is_target, model_year, NA)
    }else{
      rc_list[[rcs]][rc_list[[rcs]]$year == model_year, ]$latest_target_year <- rc_list[[rcs]][rc_list[[rcs]]$year == model_year-1, ]$latest_target_year
      if(any(rc_list[[rcs]]$year == model_year & rc_list[[rcs]]$is_target)){
        rc_list[[rcs]][rc_list[[rcs]]$year == model_year & rc_list[[rcs]]$is_target, ]$latest_target_year <- model_year}
    }
  }

  return(rc_list)               
}

#' @name get_confirmed_incidence_rate
#' @title get_confirmed_incidence_rate
#' @description to get the confirmed incidence rate in the rc_list 
#' @param rc_list
#' @param model_year
#' @param surveillance_scenario
#' @param omicron_dataset
#' @param rcs
#' @param use_mean_ir 
#' @param mean_ir_span
#' @param testing_sensitivity
#' @return 
#' @export
#' @include
get_confirmed_incidence_rate <- function(rc_list, model_year, surveillance_scenario, omicron_dataset, rcs, use_mean_ir, mean_ir_span, testing_sensitivity){
  ##### The same district/country can only use the same confirmation rate across years but different across layers 

  ### Get the confirmation rate
  rc_list[[rcs]][rc_list[[rcs]]$year == model_year, ]$confirmation_lens <- surveillance_scenario
  
  if(surveillance_scenario == "no-estimate"){ #clinical cases = confirmed cases
    no_estimate <- testing_sensitivity * 1
    rc_list[[rcs]][rc_list[[rcs]]$year == model_year, ]$confirmation_rate <- no_estimate

  }else if(surveillance_scenario == "global-estimate"){
    ## global estimate scenario 
    if(all(is.na(rc_list[[rcs]]$confirmation_rate))){
      # global_estimate <- rnorm(n = 1, mean = omicron_dataset$mean, sd = omicron_dataset$sd)
      global_estimate <- weighted.mean(rc_list$rc2[rc_list$rc2$year == model_year, ]$true_confirm_rate, 
        rc_list$rc2[rc_list$rc2$year == model_year, ]$true_incidence_rate * rc_list$rc2[rc_list$rc2$year == model_year, ]$pop_model / rc_list$rc2[rc_list$rc2$year == model_year, ]$true_confirm_rate) #just admin2 because the true confirm rate is at admin2 level, weighted on the observed cases
      global_estimate <- testing_sensitivity * global_estimate
    }else{
      global_estimate <- unique(rc_list[[rcs]]$confirmation_rate)[1]
    }
    rc_list[[rcs]][rc_list[[rcs]]$year == model_year, ]$confirmation_rate <- global_estimate
    
  }else if(surveillance_scenario == "district-estimate"){
    ## district estimate scenario 
    district_estimate <- testing_sensitivity * rc_list[[rcs]][rc_list[[rcs]]$year == model_year, ]$true_confirm_rate
    rc_list[[rcs]][rc_list[[rcs]]$year == model_year, ]$confirmation_rate <- district_estimate
  }

  ### Update the confirmed_incidence_rate
  rc_list[[rcs]][rc_list[[rcs]]$year == model_year, ]$confirmed_incidence_rate <- rc_list[[rcs]][rc_list[[rcs]]$year == model_year, ]$true_incidence_rate / rc_list[[rcs]][rc_list[[rcs]]$year == model_year, ]$true_confirm_rate * rc_list[[rcs]][rc_list[[rcs]]$year == model_year, ]$confirmation_rate
  ## get the mean_confirmed_incidence_rate
  for(i in 1:mean_ir_span){
    if(i == 1){
      mean_confirmed_incidence_rate_tmp <- rc_list[[rcs]][rc_list[[rcs]]$year == model_year, ]$confirmed_incidence_rate
    }else if(model_year - min(rc_list[[rcs]]$year) + 1 >= i){
      mean_confirmed_incidence_rate_tmp <- cbind( mean_confirmed_incidence_rate_tmp, 
                                                  rc_list[[rcs]][rc_list[[rcs]]$year == model_year - (i-1), ]$confirmed_incidence_rate )
    }else{
      mean_confirmed_incidence_rate_tmp <- cbind( mean_confirmed_incidence_rate_tmp, 
                                                  rc_list[[rcs]][rc_list[[rcs]]$year == min(rc_list[[rcs]]$year), ]$confirmed_incidence_rate )
    }
  }
  mean_confirmed_incidence_rate_tmp <- apply(mean_confirmed_incidence_rate_tmp, 1, mean)
  rc_list[[rcs]][rc_list[[rcs]]$year == model_year, ]$mean_confirmed_incidence_rate <- mean_confirmed_incidence_rate_tmp 
  rm(mean_confirmed_incidence_rate_tmp)

  return(rc_list)
}

#' @name surveillance_add_rc_new_row
#' @title surveillance_add_rc_new_row
#' @description quite similar to the load baseline incidence function, but focused on using the expected cases raster and adding new rows instead 
#' @param rc_list
#' @param ec_list
#' @param pop
#' @param model_year
#' @param sim_start_year
#' @param sim_end_year
#' @param shp1
#' @param shp2
#' @param nsamples
#' @param mean_ir_span
#' @return 
#' @export
#' @include
surveillance_add_rc_new_row <- function(rc_list, ec_list, pop, model_year, sim_start_year, sim_end_year, shp1, shp2, nsamples, mean_ir_span){

  ## Prepare
  message("Start adding new rows now. ")
  if(!is.null(ec_list$ec_rasterStack_admin1)){admin1_lambda_list <- ec_list$ec_rasterStack_admin1 / pop}else{admin1_lambda_list <- NULL}
  if(!is.null(ec_list$ec_rasterStack_admin2)){admin2_lambda_list <- ec_list$ec_rasterStack_admin2 / pop}else{admin2_lambda_list <- NULL}
  pop1 <- get_admin_population(pop, shp1) 
  pop2 <- get_admin_population(pop, shp2)
  
  # if(!is.null(ec_list$ec_rasterStack_admin1)){
  #   true_case1_all <- list()
  #   for(ly_idx in 1:nsamples){
  #     true_case1_single <- exactextractr::exact_extract(ec_list$ec_rasterStack_admin1[[ly_idx]], shp1, 'sum')
  #     true_case1_all <- append(true_case1_all, true_case1_single)
  #   }
  # }else{true_case1_all <- NULL}

  # if(!is.null(ec_list$ec_rasterStack_admin2)){
  #   true_case2_all <- list()
  #   for(ly_idx in 1:nsamples){
  #     true_case2_single <- exactextractr::exact_extract(ec_list$ec_rasterStack_admin2[[ly_idx]], shp2, 'sum')
  #     true_case2_all <- append(true_case2_all, true_case2_single)
  #   }
  # }else{true_case2_all <- NULL}
  # if(!is.null(ec_list$ec_rasterStack_admin1)){true_case1_all <- raster::extract(ec_list$ec_rasterStack_admin1, shp1, method = 'bilinear', fun = sum, na.rm = TRUE) %>% as.data.frame()}else{true_case1_all <- NULL}
  # if(!is.null(ec_list$ec_rasterStack_admin2)){true_case2_all <- raster::extract(ec_list$ec_rasterStack_admin2, shp2, method = 'bilinear', fun = sum, na.rm = TRUE) %>% as.data.frame()}else{true_case2_all <- NULL}
  # if(!is.null(ec_list$ec_rasterStack_admin1)){true_case1_all <- exactextractr::exact_extract(ec_list$ec_rasterStack_admin1, shp1, 'sum')}else{true_case1_all <- NULL}
  # if(!is.null(ec_list$ec_rasterStack_admin2)){true_case2_all <- exactextractr::exact_extract(ec_list$ec_rasterStack_admin2, shp2, 'sum')}else{true_case2_all <- NULL}

  # new try
  if(!is.null(ec_list$ec_rasterStack_admin1)){true_case1_all <- case_by_district(ec_list$ec_rasterStack_admin1, shp1)}else{true_case1_all <- NULL}
  if(!is.null(ec_list$ec_rasterStack_admin2)){true_case2_all <- case_by_district(ec_list$ec_rasterStack_admin2, shp2)}else{true_case2_all <- NULL}
  rm(ec_list)
  
  # total population
  total_pop <- sum(pop1)
  # total population test
  if(sum(pop1) != sum(pop2)){warning(paste(" *** population summarized across admin1 and admin2 districts differ by", abs(sum(pop1) - sum(pop2))))}
  

  ## loop through all the layers 
  for (layer_idx in 1:nsamples){

    if(!is.null(admin1_lambda_list)){lambda1 <- admin1_lambda_list[[layer_idx]]; lambda1[is.infinite(lambda1)] <- 0}
    if(!is.null(admin2_lambda_list)){lambda2 <- admin2_lambda_list[[layer_idx]]; lambda2[is.infinite(lambda2)] <- 0}
    if(!is.null(true_case1_all)){true_case1 <- true_case1_all[[layer_idx]]}
    if(!is.null(true_case2_all)){true_case2 <- true_case2_all[[layer_idx]]}

    ## fill in the true cases for the simulated year and end here if this is the last year of simulation 
    if(!is.null(true_case1_all)){rc_list[[layer_idx]]$rc1$true_case[rc_list[[layer_idx]]$rc1$year == model_year] <- true_case1}
    if(!is.null(true_case2_all)){rc_list[[layer_idx]]$rc2$true_case[rc_list[[layer_idx]]$rc2$year == model_year] <- true_case2}
    if(model_year >= sim_end_year){next} #************************************************************** exit when all the simulations are done
    

    if(!is.null(admin1_lambda_list)){
      ## summarize rasters to admin level 1 and copy the previous confirm rate
      incid1 <- pop_weighted_admin_mean_incid(datapath, modelpath, incidence_rate_raster = lambda1, pop_raster = pop, country, admin_shp = shp1)
      true_confirm_rate_admin1 <- head(rc_list[[layer_idx]]$rc1, nrow(shp1))$true_confirm_rate 

      # get incidence and pop table (admin 1 level)
      rc1 <- dplyr::mutate(shp1, 
                          true_incidence_rate = incid1,
                          confirmation_lens = NA,
                          confirmation_rate = NA, 
                          confirmed_incidence_rate = NA,  
                          mean_confirmed_incidence_rate = NA,  
                          pop_model = pop1,
                          pop_prop = pop1/total_pop,
                          year = model_year + 1, 
                          latest_target_year = NA,
                          is_target = NA,  # whether that place was target in this year
                          actual_prop_vaccinated = NA,
                          actual_fvp = NA, 
                          true_confirm_rate = true_confirm_rate_admin1) %>% # fully vaccinated population
        dplyr::mutate(  true_incidence_rate = ifelse(is.na(true_incidence_rate), 0, true_incidence_rate), 
                        true_case = NA ) %>% 
        dplyr::select(ISO, NAME_0, NAME_1, true_confirm_rate, true_incidence_rate, 
                      confirmation_lens, confirmation_rate, confirmed_incidence_rate, mean_confirmed_incidence_rate,  
                      pop_model, pop_prop, year, latest_target_year, is_target,
                      actual_prop_vaccinated, actual_fvp, true_case)
      
      # get the mean_true_incidence_rate
      for(i in 1:mean_ir_span){
        if(i == 1){
          mean_true_incidence_rate_tmp <- rc1$true_incidence_rate
        }else if(model_year - sim_start_year + 1 >= i){
          mean_true_incidence_rate_tmp <- cbind( mean_true_incidence_rate_tmp, rc_list[[layer_idx]]$rc1[rc_list[[layer_idx]]$rc1$year == model_year-(i-2), ]$true_incidence_rate )
        }else{
          mean_true_incidence_rate_tmp <- cbind( mean_true_incidence_rate_tmp, rc_list[[layer_idx]]$rc1[rc_list[[layer_idx]]$rc1$year == sim_start_year, ]$true_incidence_rate )
        }
      }
      mean_true_incidence_rate_tmp <- apply(mean_true_incidence_rate_tmp, 1, mean)
      rc1$mean_true_incidence_rate <- mean_true_incidence_rate_tmp 
      rm(mean_true_incidence_rate_tmp)
      rc1 %>%
        dplyr::select(ISO, NAME_0, NAME_1, true_confirm_rate, true_incidence_rate, mean_true_incidence_rate, 
                      confirmation_lens, confirmation_rate, confirmed_incidence_rate, mean_confirmed_incidence_rate,  
                      pop_model, pop_prop, year, latest_target_year, is_target,
                      actual_prop_vaccinated, actual_fvp, true_case)

    }else{rc1 <- NULL}

    if(!is.null(admin2_lambda_list)){
      ## summarize rasters to admin level 2 and copy the previous confirm rate
      incid2 <- pop_weighted_admin_mean_incid(datapath, modelpath, incidence_rate_raster = lambda2, pop_raster = pop, country, admin_shp = shp2)
      true_confirm_rate_admin2 <- head(rc_list[[layer_idx]]$rc2, nrow(shp2))$true_confirm_rate    
      
      # get incidence and pop table (admin 2 level)
      rc2 <- dplyr::mutate(shp2, 
                          true_incidence_rate = incid2,
                          confirmation_lens = NA,
                          confirmation_rate = NA, 
                          confirmed_incidence_rate = NA, 
                          mean_confirmed_incidence_rate = NA,   
                          pop_model = pop2,
                          pop_prop = pop2/total_pop,
                          year = model_year + 1, 
                          latest_target_year = NA, 
                          is_target = NA,     # whether that place was target in this year
                          actual_prop_vaccinated = NA,
                          actual_fvp = NA, 
                          true_confirm_rate = true_confirm_rate_admin2) %>% # fully vaccinated population
        dplyr::mutate(  true_incidence_rate = ifelse(is.na(true_incidence_rate), 0, true_incidence_rate), 
                        true_case = NA ) %>% 
        dplyr::select(ISO, NAME_0, NAME_1, NAME_2, true_confirm_rate, true_incidence_rate, 
                      confirmation_lens, confirmation_rate, confirmed_incidence_rate, mean_confirmed_incidence_rate,  
                      pop_model, pop_prop, year, latest_target_year, is_target,
                      actual_prop_vaccinated, actual_fvp, true_case)
      
      # get the mean_true_incidence_rate
      for(i in 1:mean_ir_span){
        if(i == 1){
          mean_true_incidence_rate_tmp <- rc2$true_incidence_rate
        }else if(model_year - sim_start_year + 1 >= i){
          mean_true_incidence_rate_tmp <- cbind( mean_true_incidence_rate_tmp, rc_list[[layer_idx]]$rc2[rc_list[[layer_idx]]$rc2$year == model_year-(i-2), ]$true_incidence_rate )
        }else{
          mean_true_incidence_rate_tmp <- cbind( mean_true_incidence_rate_tmp, rc_list[[layer_idx]]$rc2[rc_list[[layer_idx]]$rc2$year == sim_start_year, ]$true_incidence_rate )
        }
      }
      mean_true_incidence_rate_tmp <- apply(mean_true_incidence_rate_tmp, 1, mean)
      rc2$mean_true_incidence_rate <- mean_true_incidence_rate_tmp 
      rm(mean_true_incidence_rate_tmp)
      rc2 %>%
        dplyr::select(ISO, NAME_0, NAME_1, true_confirm_rate, true_incidence_rate, mean_true_incidence_rate, 
                      confirmation_lens, confirmation_rate, confirmed_incidence_rate, mean_confirmed_incidence_rate,  
                      pop_model, pop_prop, year, latest_target_year, is_target,
                      actual_prop_vaccinated, actual_fvp, true_case)

    }else{rc2 <- NULL}

    # check on population proportion 
    if ((!is.null(rc1) & !is.null(rc2)) & (sum(rc1$pop_prop)!=1 | sum(rc2$pop_prop)!=1)){
      warning(paste(" *** The population proportion calculation is incorrect for", country, ", they are", sum(rc1$pop_prop), "and", sum(rc2$pop_prop)))
    }
    
    # sometimes rc1 or rc2 has NA's for suspected_incidence due to disagreement about NA's in grid cells 
    if(!is.null(rc1)){rc_list[[layer_idx]]$rc1 <- rbind(rc_list[[layer_idx]]$rc1, rc1)}
    if(!is.null(rc2)){rc_list[[layer_idx]]$rc2 <- rbind(rc_list[[layer_idx]]$rc2, rc2)}
    rm(lambda1, lambda2, rc1, rc2, true_confirm_rate_admin1, true_confirm_rate_admin2)
  }

  rm(pop, pop1, pop2, admin1_lambda_list, admin2_lambda_list)
  gc()
  return(rc_list)
  
}
