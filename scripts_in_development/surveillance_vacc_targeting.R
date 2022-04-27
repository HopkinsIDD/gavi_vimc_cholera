# This is the draft code for VIMC model updates

# Main changes: the updated model should incorporate incidence reduction after every round of campaign
# so it is generally repetitions of:
# 0) prepare baseline incidence tables (admin1 & admin2);
# 1) select places to target and mark which areas are targeted and which are not in that year;
# 2) get the number n: this is the n-th year after the latest campaign in an area (represented by n_years_since_last_campaign column in the incidence tables);
# 3) get raster of susceptibles/incidence (considering direct and indirect effect of the latest campaign), 
#    based on the last year's incidence and n_years_since_last_campaign+1 (time elapse since the latest vacc campaign, might need this to get the residual vaccine effect)
# 4) update/calculate incidence tables based on the susceptible raster, which is the input for the next round of targeting
# ... (repeat the process every year)


#### PENDING: need to incorporate the proportion positive
# confirmation rate is a raster? or it is vector with one value for each admin ? might need to use this to calculate the true incidence (the true incidence will be compared to the threshold)

# Function to load admin0 level shp
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

# Function to load admin1 level shp
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

# function to load admin2 level shp
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



# function to load baseline incidence of a country at both admin1 and admin2 level
# return a list containing two tables(one for admin1 and one for admin2)
# each table contains admin info, baseline incidence, and pop_prop

load_baseline_incidence <- function(datapath, 
                                    country, 
                                    campaign_cov, 
                                    baseline_year, 
                                    first_vacc_year, 
                                    incidence_rate_trend, 
                                    use_country_incid_trend, 
                                    shp0, 
                                    shp1, 
                                    shp2){
  
  ## incidence data ##
  message(paste0("Loading ", datapath, "/incidence/afro_2010-2016_lambda_5k_mean.tif"))
  afr <- raster::raster(paste0(datapath, "/incidence/afro_2010-2016_lambda_5k_mean.tif"))
  country_baseline <- raster::mask(raster::crop(raster::stack(afr), shp0, snap = "out"), shp0, updatevalue = NA)
  rm(afr)
  # load population data of the first year 
  pop_baseline <- ocvImpact::create_model_pop_raster(datapath, modelpath, country, baseline_year)
  # adjust country_baseline so that it's in line with the expected cases rasters in the future (more NA cells after multiplying it with pop raster)
  pop_baseline[!is.na(raster::values(pop_baseline)), ] <- 1
  country_baseline <- country_baseline * pop_baseline


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
  
  
  ## summarize rasters to admin level 1
  incid1 <- pop_weighted_admin_mean_incid(datapath, modelpath, incidence_rate_raster = country_baseline, pop_raster = pop_baseline, country, admin_shp = shp1)
  pop1 <- get_admin_population(pop_baseline, shp1) 
  
  ## summarize rasters to admin level 2
  incid2 <- pop_weighted_admin_mean_incid(datapath, modelpath, incidence_rate_raster = country_baseline, pop_raster = pop_baseline, country, admin_shp = shp2)
  pop2 <- get_admin_population(pop_baseline, shp2)
  
  # total population
  total_pop <- sum(pop1)

  # total population test
  if(sum(pop1) != sum(pop2)){warning(paste(" *** population summarized across admin1 and admin2 districts differ by", abs(sum(pop1) - sum(pop2))))}
  
  # get incidence and pop table (admin 1 level)
  rc1 <- dplyr::mutate(shp1, 
                       suspected_incidence = incid1,
                       confirmation_lens = NA,
                       confirmation_rate = NA, 
                       confirmed_incidence = NA,  
                       pop_model = pop1,
                       pop_prop = pop1/total_pop,
                       year = baseline_year, 
                       latest_target_year = NA,
                       is_target = NA,  # whether that place was target in this year
                       actual_prop_vaccinated = NA,
                       actual_fvp = NA) %>% # fully vaccinated population
    dplyr::select(ISO, NAME_0, NAME_1, suspected_incidence, confirmation_lens, confirmation_rate, confirmed_incidence, 
                  pop_model, pop_prop, year, latest_target_year, is_target,
                  actual_prop_vaccinated, actual_fvp)
  
  # get incidence and pop table (admin 2 level)
  rc2 <- dplyr::mutate(shp2, 
                       suspected_incidence = incid2,
                       confirmation_lens = NA,
                       confirmation_rate = NA, 
                       confirmed_incidence = NA,  
                       pop_model = pop2,
                       pop_prop = pop2/total_pop,
                       year = baseline_year, 
                       latest_target_year = NA, 
                       is_target = NA,     # whether that place was target in this year
                       actual_prop_vaccinated = NA,
                       actual_fvp = NA) %>% # fully vaccinated population
    dplyr::select(ISO, NAME_0, NAME_1, NAME_2, suspected_incidence, confirmation_lens, confirmation_rate, confirmed_incidence, 
                  pop_model, pop_prop, year, latest_target_year, is_target,
                  actual_prop_vaccinated, actual_fvp)
  
  # check on population proportion 
  if (sum(rc1$pop_prop)!=1 | sum(rc2$pop_prop)!=1){
   warning(paste(" *** The population proportion calculation is incorrect for", country, ", they are", sum(rc1$pop_prop), "and", sum(rc2$pop_prop)))
  }
  
  rc_list <- list("rc1" = rc1, "rc2" = rc2)
  rm(pop1, pop2)
  gc()
  
  return(rc_list)
  
}



# Function to update target indicator (is_target column) and vaccine coverage (actual_prop_vaccinated column)
# is_target: to mark whether a place in a certain year is targeted or not 
# actual_prop_vaccinated: if that place is vaccinated, then this is campaign_cov = 0.8, if it that place is not vaccinated, then it is 0.
update_targets_list <- function(datapath, modelpath, country, scenario, 
                                rc_list, model_year, #here the baseline year is the first vacc year, but it's not likely to be used 
                                campaign_cov, 
                                threshold, vac_unconstrained, 
                                surveillance_scenario, 
                                vac_interval, #this is country level 
                                vac_start_year, vac_end_year, 
                                num_skip_years #this is district level 
                                ){
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
      ( !all(is.na(rc_list$rc1$latest_target_year)) & (model_year - max(rc_list$rc1$latest_target_year, na.rm = TRUE) < vac_interval 
                                                        | model_year - max(rc_list$rc2$latest_target_year, na.rm = TRUE) < vac_interval) ) | 
      (scenario == "no-vaccination") |
      (model_year > vac_end_year) |
      (model_year < vac_start_year) ){
    
    message(paste("The vaccination campaign is skipped for year", model_year, "for the whole country", country))
    if( max(rc_list$rc1$latest_target_year) != max(rc_list$rc2$latest_target_year) ){
      stop("admin1 and admin2 do not agree on the latest target year in the rc list table, please check. ")
    }else{
      rc_list <- update_novacc_year(rc_list, model_year)
    }

  }else{
    rc_list <- update_vacc_year(datapath, modelpath, country, rc_list, model_year, 
                                campaign_cov, threshold, surveillance_scenario, 
                                vac_start_year, vac_end_year, num_skip_years) 
  }
  
  return(rc_list) # updated with latest targets of this year
}



update_novacc_year <- function(rc_list, model_year){
  # Update the latest target year, if all NA's from last year or this is the first simulation year, no need to update
  if(model_year > min(rc_list$rc1$year)){ 
    if(!all(is.na(rc_list$rc1[rc_list$rc1$year == model_year-1, ]$latest_target_year))){
      rc_list$rc1[rc_list$rc1$year == model_year, ]$latest_target_year <- rc_list$rc1[rc_list$rc1$year == model_year-1, ]$latest_target_year
    }
    if(!all(is.na(rc_list$rc2[rc_list$rc2$year == model_year-1, ]$latest_target_year))){
      rc_list$rc2[rc_list$rc2$year == model_year, ]$latest_target_year <- rc_list$rc2[rc_list$rc2$year == model_year-1, ]$latest_target_year
    }
  }
  # Update the other variables 
  rc_list$rc1[rc_list$rc1$year == model_year, c("is_target", "actual_prop_vaccinated", "actual_fvp")] <- 0 
  rc_list$rc2[rc_list$rc2$year == model_year, c("is_target", "actual_prop_vaccinated", "actual_fvp")] <- 0
  return(rc_list)
  ### this is where we might need to include the waning effects of vaccines, but we're not doing it for the time being
}



update_vacc_year<- function(datapath, modelpath, country, rc_list, model_year, 
                            campaign_cov, threshold, surveillance_scenario, 
                            vac_start_year, vac_end_year, num_skip_years){
  ### Put on the confirmation lens 
  ## Read in
  message(paste0("Loading ", datapath, "/confirmation_rate/something.csv"))
  omicron_dataset <- readr::read_csv(paste0(datapath, "/confirmation_rate/something.csv")) #the true confirmation rate 
  ## Assign
  rc_list <- get_confirmed_incidence_rate(rc_list, model_year, surveillance_scenario, omicron_dataset)

  ### Update the rest of the variables 
  ## Update the target and vaccinated pop
  #if there is no previous year
  if(min(rc_list$rc1$year) == model_year){
    rc_list$rc1[rc_list$rc1$year == model_year, ]$is_target <- ifelse(rc_list$rc1[rc_list$rc1$year == model_year, ]$confirmed_incidence >= threshold, 1, 0)
    rc_list$rc2[rc_list$rc2$year == model_year, ]$is_target <- ifelse(rc_list$rc2[rc_list$rc2$year == model_year, ]$confirmed_incidence >= threshold, 1, 0)
  #if there is previous year 
  }else{
    rc_list$rc1[rc_list$rc1$year == model_year, ]$is_target <- ifelse(  (!is.na(rc_list$rc1[rc_list$rc1$year == model_year-1, ]$latest_target_year)
                                                                          & model_year - rc_list$rc1[rc_list$rc1$year == model_year-1, ]$latest_target_year >= num_skip_years
                                                                          & rc_list$rc1[rc_list$rc1$year == model_year, ]$confirmed_incidence >= threshold)
                                                                      | (is.na(rc_list$rc1[rc_list$rc1$year == model_year-1, ]$latest_target_year)
                                                                          & rc_list$rc1[rc_list$rc1$year == model_year, ]$confirmed_incidence >= threshold), 1, 0)
    rc_list$rc2[rc_list$rc2$year == model_year, ]$is_target <- ifelse(  (!is.na(rc_list$rc2[rc_list$rc2$year == model_year-1, ]$latest_target_year)
                                                                          & model_year - rc_list$rc2[rc_list$rc2$year == model_year-1, ]$latest_target_year >= num_skip_years
                                                                          & rc_list$rc2[rc_list$rc2$year == model_year, ]$confirmed_incidence >= threshold)
                                                                      | (is.na(rc_list$rc2[rc_list$rc2$year == model_year-1, ]$latest_target_year)
                                                                          & rc_list$rc2[rc_list$rc2$year == model_year, ]$confirmed_incidence >= threshold), 1, 0)
  }
  rc_list$rc1[rc_list$rc1$year == model_year, ]$actual_prop_vaccinated <- rc_list$rc1[rc_list$rc1$year == model_year, ]$is_target * campaign_cov
  rc_list$rc1[rc_list$rc1$year == model_year, ]$actual_fvp <- rc_list$rc1[rc_list$rc1$year == model_year, ]$actual_prop_vaccinated * rc_list$rc1[rc_list$rc1$year == model_year, ]$pop_model
  rc_list$rc2[rc_list$rc2$year == model_year, ]$actual_prop_vaccinated <- rc_list$rc2[rc_list$rc2$year == model_year, ]$is_target * campaign_cov
  rc_list$rc2[rc_list$rc2$year == model_year, ]$actual_fvp <- rc_list$rc2[rc_list$rc2$year == model_year, ]$actual_prop_vaccinated * rc_list$rc2[rc_list$rc2$year == model_year, ]$pop_model
  
  ## Update the latest target year
  if(min(rc_list$rc1$year) == model_year){
    rc_list$rc1[rc_list$rc1$year == model_year, ]$latest_target_year <- ifelse(rc_list$rc1[rc_list$rc1$year == model_year, ]$is_target, model_year, NA)
    rc_list$rc2[rc_list$rc2$year == model_year, ]$latest_target_year <- ifelse(rc_list$rc2[rc_list$rc2$year == model_year, ]$is_target, model_year, NA)
  }else{
    rc_list$rc1[rc_list$rc1$year == model_year, ]$latest_target_year <- rc_list$rc1[rc_list$rc1$year == model_year-1, ]$latest_target_year
    rc_list$rc1[rc_list$rc1$year == model_year & rc_list$rc1$is_target, ]$latest_target_year <- model_year
    rc_list$rc2[rc_list$rc2$year == model_year, ]$latest_target_year <- rc_list$rc2[rc_list$rc2$year == model_year-1, ]$latest_target_year
    rc_list$rc2[rc_list$rc2$year == model_year & rc_list$rc2$is_target, ]$latest_target_year <- model_year
  }

  return(rc_list)               
}



get_confirmed_incidence_rate <- function(rc_list, model_year, surveillance_scenario, omicron_dataset){
  ### Get the confirmation rate
  rc_list$rc1[rc_list$rc1$year == model_year, ]$confirmation_lens <- surveillance_scenario
  rc_list$rc2[rc_list$rc2$year == model_year, ]$confirmation_lens <- surveillance_scenario

  if(surveillance_scenario == "no-estimate"){
    rc_list$rc1[rc_list$rc1$year == model_year, ]$confirmation_rate <- 1
    rc_list$rc2[rc_list$rc2$year == model_year, ]$confirmation_rate <- 1
  }else if(surveillance_scenario == "global-estimate"){
    ### determine how to use the omicron dataset
    global_estimate <- truncnorm::rtruncnorm(n = 1, a = 0, b = 1, mean = 0.43, sd = 0.5) #tmp, for now
    rc_list$rc1[rc_list$rc1$year == model_year, ]$confirmation_rate <- global_estimate
    rc_list$rc2[rc_list$rc2$year == model_year, ]$confirmation_rate <- global_estimate
  }else if(surveillance_scenario == "district-estimate"){
    ### determine how to use the omicron dataset
    district_estimate_rc1 <- truncnorm::rtruncnorm(n = nrow(rc_list$rc1[rc_list$rc1$year == model_year, ]), a = 0, b = 1, mean = 0.43, sd = 0.5)
    rc_list$rc1[rc_list$rc1$year == model_year, ]$confirmation_rate <- district_estimate_rc1
    district_estimate_rc2 <- truncnorm::rtruncnorm(n = nrow(rc_list$rc2[rc_list$rc2$year == model_year, ]), a = 0, b = 1, mean = 0.43, sd = 0.5)
    rc_list$rc2[rc_list$rc2$year == model_year, ]$confirmation_rate <- district_estimate_rc2
  }

  ### Update the confirmed rate 
  rc_list$rc1[rc_list$rc1$year == model_year, ]$confirmed_incidence <- rc_list$rc1[rc_list$rc1$year == model_year, ]$confirmation_rate * rc_list$rc1[rc_list$rc1$year == model_year, ]$suspected_incidence
  rc_list$rc2[rc_list$rc2$year == model_year, ]$confirmed_incidence <- rc_list$rc2[rc_list$rc2$year == model_year, ]$confirmation_rate * rc_list$rc2[rc_list$rc2$year == model_year, ]$suspected_incidence
  
  return(rc_list)
}



surveillance_add_rc_new_row <- function(rc_list, ec_list, pop, model_year, sim_start_year, sim_end_year, shp1, shp2){

  ## load newly updated incidence data 
  year_idx <- match(model_year, sim_start_year:sim_end_year)
  lambda1 <- ec_list$ec_rasterStack_admin1[[year_idx]] / pop
  lambda2 <- ec_list$ec_rasterStack_admin2[[year_idx]] / pop
  
  ## summarize rasters to admin level 1
  incid1 <- pop_weighted_admin_mean_incid(datapath, modelpath, incidence_rate_raster = lambda1, pop_raster = pop, country, admin_shp = shp1)
  pop1 <- get_admin_population(pop, shp1) 
  
  ## summarize rasters to admin level 2
  incid2 <- pop_weighted_admin_mean_incid(datapath, modelpath, incidence_rate_raster = lambda2, pop_raster = pop, country, admin_shp = shp2)
  pop2 <- get_admin_population(pop, shp2)
  
  # total population
  total_pop <- sum(pop1)

  # total population test
  if(sum(pop1) != sum(pop2)){warning(paste(" *** population summarized across admin1 and admin2 districts differ by", abs(sum(pop1) - sum(pop2))))}
  
  # get incidence and pop table (admin 1 level)
  rc1 <- dplyr::mutate(shp1, 
                       suspected_incidence = incid1,
                       confirmation_lens = NA,
                       confirmation_rate = NA, 
                       confirmed_incidence = NA,  
                       pop_model = pop1,
                       pop_prop = pop1/total_pop,
                       year = model_year + 1, 
                       latest_target_year = NA,
                       is_target = NA,  # whether that place was target in this year
                       actual_prop_vaccinated = NA,
                       actual_fvp = NA) %>% # fully vaccinated population
    dplyr::select(ISO, NAME_0, NAME_1, suspected_incidence, confirmation_lens, confirmation_rate, confirmed_incidence, 
                  pop_model, pop_prop, year, latest_target_year, is_target,
                  actual_prop_vaccinated, actual_fvp)
  
  # get incidence and pop table (admin 2 level)
  rc2 <- dplyr::mutate(shp2, 
                       suspected_incidence = incid2,
                       confirmation_lens = NA,
                       confirmation_rate = NA, 
                       confirmed_incidence = NA,  
                       pop_model = pop2,
                       pop_prop = pop2/total_pop,
                       year = model_year + 1, 
                       latest_target_year = NA, 
                       is_target = NA,     # whether that place was target in this year
                       actual_prop_vaccinated = NA,
                       actual_fvp = NA) %>% # fully vaccinated population
    dplyr::select(ISO, NAME_0, NAME_1, NAME_2, suspected_incidence, confirmation_lens, confirmation_rate, confirmed_incidence, 
                  pop_model, pop_prop, year, latest_target_year, is_target,
                  actual_prop_vaccinated, actual_fvp)
  
  # check on population proportion 
  if (sum(rc1$pop_prop)!=1 | sum(rc2$pop_prop)!=1){
   warning(paste(" *** The population proportion calculation is incorrect for", country, ", they are", sum(rc1$pop_prop), "and", sum(rc2$pop_prop)))
  }
  
  # sometimes rc1 or rc2 has NA's for suspected_incidence due to disagreement about NA's in grid cells 
  rc_list$rc1 <- rbind(rc_list$rc1, rc1)
  rc_list$rc2 <- rbind(rc_list$rc2, rc2)
  rm(pop1, pop2)
  gc()
  
  return(rc_list)
  
}
