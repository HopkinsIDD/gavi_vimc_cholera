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

load_baseline_incidence <- function(datapath, country, campaign_cov = 0.8, baseline_year = 2018,
                                    shp1, shp2){
  
  ## incidence data ##
  message(paste0("Loading ", datapath, "/incidence/afro_2010-2016_lambda_5k_mean.tif"))
  afr <- raster::raster(paste0(datapath, "/incidence/afro_2010-2016_lambda_5k_mean.tif"))
  
  # load population data of the first year 
  pop_baseline <- create_model_pop_raster(datapath, modelpath, country, baseline_year)
  
  ## admin1 unit shapefile ##
  # shp1 <- load_shp1_by_country(datapath, country) # can also load the shpfiles outside this function
  # shp2 <- load_shp2_by_country(datapath, country)
  
  ## summarize rasters to admin level 1
  incid1 <- exactextractr::exact_extract(afr, shp1, 'mean') ## could add population weight here for better incidence estimate but need to project population to the incidence grid
  pop1 <- get_admin_population(pop_baseline, shp1)
  
  ## summarize rasters to admin level 2
  incid2 <- exactextractr::exact_extract(afr, shp2, 'mean') ## could add population weight here for better incidence estimate but need to project population to the incidence grid
  pop2 <- get_admin_population(pop_baseline, shp2)
  
  # total population
  total_pop <- sum(pop1)
  
  # get incidence and pop table (admin 1 level)
  rc1 <- dplyr::mutate(shp1, 
                       incidence = incid1,
                       pop_model = pop1,
                       pop_prop = pop1/total_pop,
                       target_year = baseline_year, # baseline incidence is inc of year 1
                       is_target = NA,  # whether that place was target in this year
                       actual_prop_vaccinated = NA,
                       campaign_cov = campaign_cov,
                       actual_fvp = NA) %>% # fully vaccinated population
    # sf::st_drop_geometry() %>%
    dplyr::select(ISO, NAME_0, NAME_1, incidence, pop_model, pop_prop, target_year, is_target, 
                  actual_prop_vaccinated, actual_fvp)
    # tibble::as_tibble()
  
  # get incidence and pop table (admin 2 level)
  rc2 <- dplyr::mutate(shp2, 
                       incidence = incid2,
                       pop_model = pop2,
                       pop_prop = pop2/total_pop,
                       target_year = baseline_year,    # baseline incidence is inc of year 1
                       is_target = NA,     # whether that place was target in this year
                       actual_prop_vaccinated = NA,
                       campaign_cov = campaign_cov,
                       actual_fvp = NA) %>% # fully vaccinated population
    # sf::st_drop_geometry() %>%
    dplyr::select(ISO, NAME_0, NAME_1, NAME_2, incidence, pop_model, pop_prop, target_year, is_target,
                  actual_prop_vaccinated, actual_fvp)
    # tibble::as_tibble()
  
  # this if condition is always true, skip this for now
  # if (sum(rc1$pop_prop)!=1 | sum(rc2$pop_prop)!=1){
  #  stop(paste("The population proportion calculation is incorrect for", country))
  # }
  
  rc_list <- list("rc1" = rc1, "rc2" = rc2)
  rm(afr, pop1, pop2)
  gc()
  
  return(rc_list)
  
}


# Function to update target indicator (is_target column) and vaccine coverage (actual_prop_vaccinated column)
# is_target: to mark whether a place in a certain year is targeted or not 
# actual_prop_vaccinated: if that place is vaccinated, then this is campaign_cov = 0.8, if it that place is not vaccinated, then it is 0.
update_targets_list <- function(rc_list, threshold = 1/10000, baseline_year = baseline_year, campaign_cov = campaign_cov){
  
  this_year <- max(rc_list[[1]]$target_year) #it could be directly from the argument 
  rc1_this_year <- rc_list[[1]] %>% filter(target_year == this_year)
  rc2_this_year <- rc_list[[2]] %>% filter(target_year == this_year)

  if(this_year == baseline_year){
    
    rc1_this_year <- rc1_this_year %>% mutate(is_target = ifelse(incidence > threshold, 1, 0))
    rc2_this_year <- rc2_this_year %>% mutate(is_target = ifelse(incidence > threshold, 1, 0))
  
  }

  # if this_year < baseline_year+4, then no need to check whether that place had campaign in past three years, just need to direct select the targets
  if(this_year < baseline_year + 4 & this_year > baseline_year){
    
    # check if this place has been vaccinated in the paste three years (a district cannot be vaccinated every 3 years)
    recent_targets_admin1 <- (rc_list[[1]] %>% filter(target_year < this_year & is_target == 1))$NAME_1
    recent_targets_admin2 <- (rc_list[[2]] %>% filter(target_year < this_year & is_target == 1))$NAME_2
   
    potential_targets_admin1 <- (rc1_this_year %>% filter(incidence > threshold))$NAME_1
    potential_targets_admin2 <- (rc2_this_year %>% filter(incidence > threshold))$NAME_2

    rc1_this_year <- rc1_this_year %>%
      mutate(is_target = ifelse(NAME_1 %in% potential_targets_admin1 & NAME_1 %in% recent_targets_admin1 == FALSE, 1, 0))
    rc2_this_year <- rc2_this_year %>%
      mutate(is_target = ifelse(NAME_2 %in% potential_targets_admin2 & NAME_2 %in% recent_targets_admin2 == FALSE, 1, 0))
    
  }
  
  if(this_year >= baseline_year + 4){
    
    # check if this place has been vaccinated in the paste three years (a district cannot be vaccinated every 3 years)
    recent_targets_admin1 <- (rc_list[[1]] %>% filter(target_year < this_year & target_year >= this_year - 3 & is_target == 1))$NAME_1
    recent_targets_admin2 <- (rc_list[[2]] %>% filter(target_year < this_year & target_year >= this_year - 3 & is_target == 1))$NAME_2
   
    potential_targets_admin1 <- (rc1_this_year %>% filter(incidence > threshold))$NAME_1
    potential_targets_admin2 <- (rc2_this_year %>% filter(incidence > threshold))$NAME_2

    rc1_this_year <- rc1_this_year %>%
      mutate(is_target = ifelse((NAME_1 %in% potential_targets_admin1) & (NAME_1 %in% recent_targets_admin1 == FALSE), 1, 0))
    rc2_this_year <- rc2_this_year %>%
      mutate(is_target = ifelse((NAME_2 %in% potential_targets_admin2) & (NAME_2 %in% recent_targets_admin2 == FALSE), 1, 0))
    
  }
  
  rc1_this_year <- rc1_this_year %>%
    mutate(actual_prop_vaccinated = ifelse(is_target == 1, campaign_cov, 0))
  rc2_this_year <- rc2_this_year %>%
    mutate(actual_prop_vaccinated = ifelse(is_target == 1, campaign_cov, 0))
 

  # update target tables
  rc_list[[1]] <- rc_list[[1]] %>% filter(target_year < this_year) %>% rbind(rc1_this_year)
  rc_list[[2]] <- rc_list[[2]] %>% filter(target_year < this_year) %>% rbind(rc2_this_year)
  
  # update actual_fvp
  rc_list[[1]] <- rc_list[[1]] %>% mutate(actual_fvp = pop_model * actual_prop_vaccinated)
  rc_list[[2]] <- rc_list[[2]] %>% mutate(actual_fvp = pop_model * actual_prop_vaccinated)
  
  return(rc_list) # updated with latest targets of this year
  
}


# Function to append the new template for a new year's campaign
# three parts:
# 1) append rows and fill in with admin names
# 2) fill in target_year column (which is latest_campaign_year + 1)
# 3) fill in pop_model column (based on the raster of that year)
new_campaign_preparation <- function(rc_list, shp1, shp2){ # raster layer for the new campaign year
  
  latest_campaign_year <- max(rc_list[[1]]$target_year)
  
  message(paste0("Appending new rows in rc_list tables in preparation of next year's (", latest_campaign_year+1, ") campaign."))
  new_pop <- create_model_pop_raster(datapath, modelpath, country, latest_campaign_year + 1)
  
  # prepare empty columns
  df_new_campaign_admin1 <- rc_list[[1]] %>% 
    filter(target_year == latest_campaign_year) %>%
    mutate(incidence=NA, pop_prop=NA, pop_model=NA, target_year=latest_campaign_year+1, is_target=NA, actual_prop_vaccinated=NA, actual_fvp=NA)
  df_new_campaign_admin2 <- rc_list[[2]] %>% 
    filter(target_year == latest_campaign_year) %>%
    mutate(incidence=NA, pop_prop=NA, pop_model=NA, target_year=latest_campaign_year+1, is_target=NA, actual_prop_vaccinated=NA, actual_fvp=NA)
  
  # fill in pop_model and pop_prop columns
  # pop <- create_model_pop_raster(datapath, modelpath, country, latest_campaign_year + 1) # might need to pre-run function to get this outside this function, because in the functions creating population rasters (in order to update incidence), this will be reused
  
  df_new_campaign_admin1$pop_model <- exactextractr::exact_extract(new_pop, shp1, 'sum') # might need to pre-run functions about to get shp1 and shp2
  df_new_campaign_admin2$pop_model <- exactextractr::exact_extract(new_pop, shp2, 'sum')
  
  total_pop <- sum(df_new_campaign_admin1$pop_model)
  df_new_campaign_admin1$pop_prop <- df_new_campaign_admin1$pop_model / total_pop
  df_new_campaign_admin2$pop_prop <- df_new_campaign_admin2$pop_model / total_pop
  
  # append rows of new campaign to the working tables 
  rc_list[[1]] <- rbind(rc_list[[1]], df_new_campaign_admin1)
  rc_list[[2]] <- rbind(rc_list[[2]], df_new_campaign_admin2)
  
  return(rc_list)
  
}
