
#' @name get_admin_population
#' @title get_admin_population
#' @description Wrapper for exact_extract with a given population raster and shapefile to get admin-level worldpop population summaries for a given country. Note that the default settings use na.rm = TRUE
#' @param pop population raster
#' @param shp admin-level sf object
#' @return vector of populations by admin2 in order of admin2 names in the shapefile
#' @export
get_admin_population <- function(pop, shp){
  ## TODO: add test for GID2 and column names
  rc <- exactextractr::exact_extract(pop, shp, 'sum')
  return(rc)
}


#' @name load_targets_by_country
#' @title load_targets_by_country
#' @description Load possible vaccination targets by country, including incidence and population proportion information
#' @param datapath path to data
#' @param modelpath modelpath
#' @param country country code
#' @importFrom magrittr %>%
#' @return dataframe with incidence and population by admin unit for a single country
#' @export
#' @include load_worldpop_by_country.R load_shapefile_by_country.R utils_targeting.R
load_targets_by_country <- function(datapath, modelpath, country){
  ### prepare for the full model run
  group_id <- 'JHU-Lee'
  SplittedString = strsplit(modelpath, '/')[[1]]
  touchstone = SplittedString[length(SplittedString)] #generated from the modelpath

  ExpectationsIDList <- c()
  for (teams in montagu::montagu_expectations(group_id, touchstone)$description){
    if (group_id %in% strsplit(teams, ":")[[1]]){
      ExpectationsIDList <- c(ExpectationsIDList, montagu::montagu_expectations(group_id, touchstone)$id[match(teams, montagu::montagu_expectations(group_id, touchstone)$description)])
    }
  }
  if (length(ExpectationsIDList) == 1){
    CountriesForSim <- montagu::montagu_expectation_countries(group_id, touchstone, ExpectationsIDList)
  }else{
    message('There are multiple expectations for the current touchstone and gourp id being used, the Montagu API cannot return a single country list. ')
    message('The error is within load_targets_by_country function. ')
  }

  if (country %in% CountriesForSim$id){

    ## incidence data ##

    ##calam added 12/28/23 to use new raster for the 2023 touchstone
    runname <- config$runname
    if(runname == "202310gavi-4"){
      message(paste0("Loading ", datapath, "/incidence/afro_2016-2020_lambda_5k_mean.tif"))
      afr <- raster::raster(paste0(datapath, "/incidence/afro_2016-2020_lambda_5k_mean.tif"))
    } else {
      message(paste0("Loading ", datapath, "/incidence/afro_2010-2016_lambda_5k_mean.tif"))
      afr <- raster::raster(paste0(datapath, "/incidence/afro_2010-2016_lambda_5k_mean.tif"))
    }
    

    ##end addition

    ## WorldPop population data ##
    pop <- load_worldpop_by_country(datapath, country)
    ## admin unit shapefile ##
    
    ##The VIMC Core Model
    if (as.logical(config$custom$use_custom_shapefile) == FALSE){ 
      message("load vaccine targets using the GADM admin 2 shapefile")
      shp <- load_shapefile_by_country(datapath, country)
    } else { ## The DRC Case Study, which uses a custom shapefile for health zones
      message("load vaccine targets using the custom health zone shapefile")
      shp <- load_custom_shapefile_by_country(admin0 = FALSE)
      sf::st_crs(shp) <- 4326 ## for some reason crs needs to be re-set after loading the custom shapefile (to investigate)
    }

    ## summarize rasters to admin level (BGD, non-raster, and african raster countries)
    if (country == "BGD"){
      bgd <- raster::raster(paste0(datapath, "/incidence/BGD_incid_5k_100.tif"))
      incid2 <- exactextractr::exact_extract(bgd, shp, 'mean') ## to use BGD incidence raster for targeting
    } else if (country %in% c("AFG", "HTI", "IRN", "IRQ", "NPL", "PAK", "PHL", "THA", "YEM", "IND")) { #for non-raster countries
       number_samples <- config$incid$num_samples
       file_path <- paste0(datapath, "/incidence/",country, "_incid_5k_", number_samples, ".tif")
       rm(number_samples) ##to remove global variable
       message(paste0("Loading ", file_path))
       ##check that the incidence rate raster is cropped and exists 
       if (file.exists(file_path)){
         non_rast <- raster::raster(file_path)
         incid2 <- exactextractr::exact_extract(non_rast, shp, 'mean') ## to use incidence raster for targeting for non-raster countries
       } else {
          stop(paste("Run the run_country_incid_crop script first for ", country))
       }
    } else {
      incid2 <- exactextractr::exact_extract(afr, shp, 'mean') ## could add population weight here for better incidence estimate but need to project population to the incidence grid
    }
    
    
    ## 30 Apr 2024: test for DRC Case study: filter out health zones with NA incidence
    
    ## if(any(is.na(incid2))){
      ##print(" NA values found in loaded targets' incidence, removing units with NA incidence ")
      ##shp <- shp[-which(is.na(incid2)),] ## remove rows with NA incidence (units that fall outside the incidence raster)
      ##message(paste0(" removed units with NA incidence from shapefile: ", shp$NAME_2[which(is.na(incid2))]))
      ##incid2 <- na.omit(incid2) ## remove NAs from health zone/admin2 incidence
    ##}
    
    ## end test
    
    pop2 <- get_admin_population(pop, shp)
    total_pop <- sum(pop2)
    
    ## 30 Apr 2024 debugging - check number of rows for population per admin unit
    message(paste0("the population per health zone table has ", length(pop2), " elements"))

    ### do a little thing to the dataframe -- 7/2021
    ## This is only necessary when using the GADM admin 2 shapefile, since the custom DRC shapefile already has the required columns from shp_sp
    
    if(as.logical(config$custom$use_custom_shapefile) == FALSE){
      shp <- shp %>%
        dplyr::mutate(genID = paste0(NAME_0, '-', NAME_1, '-', NAME_2))
      shp_sp <- GADMTools::gadm_sp_loadCountries(c(country), level = 2, basefile = file.path(datapath, "shapefiles/"))$spdf
      shp_sp$genID <- paste0(shp_sp$NAME_0, '-', shp_sp$NAME_1, '-', shp_sp$NAME_2)
      shp <- merge(shp, shp_sp, id = 'genID')      
    }

    rc <- dplyr::mutate(shp,
                        incidence = incid2,
                        pop_wp = pop2,
                        pop_prop = pop2/total_pop) %>%
      sf::st_drop_geometry() %>%
      dplyr::select(GID_0, GID_2, NAME_1, NAME_2, incidence, pop_prop) %>%
      tibble::as_tibble() 
    

    if (sum(rc$pop_prop)!=1){
      stop(paste("The population proportion calculation is incorrect for", country))
    } ## DEBUG comment this section

    rm(afr,pop,shp)
    gc()

  } else{

    stop(paste(country, "vaccine targeting is not yet supported."))
  }

  return(rc)
}


#' @name assign_vaccine_targets
#' @title assign_vaccine_targets
#' @description Assign vaccine to a country and year
#' @param datapath path to external data
#' @param modelpath path to model input (montagu) data
#' @param country country code
#' @param scenario Unique string that identifies the coverage scenario name
#' @param cache montagu cache
#' @param targeting_strat Unique string that identifies how targets should be ranked (default = incidence)
#' @param campaign_cov Proportion of the target population vaccinated by the campaign (default = 0.8)
#' @param num_skip_years Skip targeting in locations vaccinated within the last X years (default = 0 --> no skipping)
#' @importFrom magrittr %>%
#' @return dataframe
#' @export
#' @include utils_targeting.R utils_montagu.R
assign_vaccine_targets <- function(datapath, modelpath, country, scenario, cache, targeting_strat = "incidence", campaign_cov = 0.8, num_skip_years = 0){

  message(paste("Now assigning vaccine by incidence:", country, scenario))
  ptargets <- load_targets_by_country(datapath, modelpath, country)
  ###########add a little check point for the situation when the coverage data exists but is just 0
  if (as.logical(config$custom$use_montagu_coverage) == TRUE){  ##the VIMC Core Model
    coverage <- import_coverage_scenario(modelpath, country, scenario, cache, filter0 = FALSE, redownload = FALSE)
  } else {  ## the DRC Case Study
    coverage <- import_coverage_scenario_custom(datapath, country, scenario, cache, filter0 = FALSE)
  }

  coverage_as_all_0_for_campaign <- (sum(coverage$OCV1) == 0) & (sum(coverage$OCV2) == 0)

  if (!coverage_as_all_0_for_campaign){
    ## this seems to be a catch-all filtering step, in case montagu starts including no-vaccination years in its coverage sheets
    
    if (as.logical(config$custom$use_montagu_coverage) == TRUE){  ##the VIMC Core Model
      coverage <- import_coverage_scenario(modelpath, country, scenario, cache, filter0 = FALSE, redownload = FALSE)
    } else {  ## the DRC Case Study
      coverage <- import_coverage_scenario_custom(datapath, country, scenario, cache, filter0 = FALSE)
    }

  }

  if (is.null(coverage)){

    message("No vaccine was assigned for the no-vaccination scenario. Returning NULL value.")
    ftargets_flat <- NULL

  } else{

    ## Perform checks on the coverage scenario
    if (!all(coverage$gender == "both") |
        !all(coverage$age_range_verbatim == "default age and gender" | coverage$age_range_verbatim == ">1y" | coverage$age_range_verbatim == "default age groups" | coverage$age_range_verbatim == "<NA>" | coverage$age_range_verbatim == "1-100")|
        !all(coverage$activity_type == "campaign")
    ){
      stop(paste("Vaccine assignment is not supported for this coverage scenario. Check the gender, age_range_verbatim, and activity_type columns in the", scenario, "coverage sheet."))
    }

    ftargets <- vector(mode = "list", length = length(coverage$year))
    print(coverage$year)
    for (i in 1:length(coverage$year)){

      cov_year <- coverage[i,]
      
      ## Major modification 1 May 2024 for DRC case study - testing
      
      if(as.logical(config$custom$use_montagu_coverage) == TRUE){
        
        ## montagu coverage tables have a target population that represents the country level population
        
        goal_target_pop <- cov_year$target ## target population for the vaccination campaign @ country level
        
      } else {
        
        ## if we are not using montagu coverage, target population needs to be made equal to country population
        ## otherwise, as in the DRC case study, since the target population in the coverage table represents a small
        ## subset of the country population, we get low numbers of possible vaccinated people in the calculation in lines 253-254
        
        ## get country population from montagu excluding people under 1
        
        goal_target_pop <- import_int_country_population(modelpath, country, cache, redownload = FALSE) %>%
          dplyr::filter(country_code == !!country) %>%
          dplyr::filter(year == coverage$year[i]) %>%
          dplyr::filter(!age_from %in% c(0)) %>%   ##this filters out the people aged 0
          dplyr::rename(pop_age = value) %>%
          dplyr::select(country_code, year, age_from, age_to, pop_age) %>%
          dplyr::summarise(pop_tot = sum(pop_age))
        
        ##debugging
        str(goal_target_pop)
        
        ## keep only the numeric value
        goal_target_pop <- as.numeric(goal_target_pop[1,1])
        
        message("turned goal_target_pop to numeric")
        
        ##debugging
        str(goal_target_pop)
        
           
        
      }

      ## end major modification
      
      ##goal fvps with one dose and two doses of the vaccine
      goal_ocv1_fvp <- cov_year$fvp_ocv1
      goal_ocv2_fvp <- cov_year$fvp_ocv2
      prop_ocv1 <- cov_year$prop_ocv1 ## proportion of vaccinees that receive 1 dose

      ## Skip locations vaccinated within the previous `num_skip_years`
      skip_years <- cov_year$year - (1:num_skip_years)
      if(num_skip_years>0 & length(skip_years)>0 & any(skip_years %in% coverage$year)){

        message(paste("Skip locations vaccinated in", paste(skip_years, collapse = ",")))
        incl_skip_years <- which(names(ftargets) %in% skip_years)
        recent_targets <- unique(as.vector(sapply(incl_skip_years, function(x){
          ftargets[[x]]$GID_2
        })))

        new_ptargets <- dplyr::filter(ptargets, !(GID_2 %in% recent_targets))

      } else{
        message(paste("Skipping was disabled (num_skip_years=0) or there were no recent campaigns within", num_skip_years, "years of", cov_year$year))
        new_ptargets <- ptargets
      }

      ## full district targeting
      ptargets_avail <- run_targeting_strategy(new_ptargets, targeting_strat) %>%
        dplyr::mutate(
          id = seq_along(GID_2),
          possible_fvp_ocv1 = round(goal_target_pop * pop_prop * campaign_cov * prop_ocv1, 0),
          possible_fvp_ocv2 = round(goal_target_pop * pop_prop * campaign_cov * (1-prop_ocv1), 0),
          possible_fvp_cumsum_ocv1 = cumsum(possible_fvp_ocv1),
          possible_fvp_cumsum_ocv2 = cumsum(possible_fvp_ocv2),
          actual_ocv1_fvp = ifelse(possible_fvp_cumsum_ocv1 <= goal_ocv1_fvp, possible_fvp_ocv1, 0),
          actual_ocv2_fvp = ifelse(possible_fvp_cumsum_ocv2 <= goal_ocv2_fvp, possible_fvp_ocv2, 0),
          fullDistrict_ocv1_target = ifelse(actual_ocv1_fvp > 0, TRUE, FALSE),
          fullDistrict_ocv2_target = ifelse(actual_ocv2_fvp > 0, TRUE, FALSE))

      ## partial district targeting
      fullDistrict_vaccinated_ocv1_people <- dplyr::filter(ptargets_avail, fullDistrict_ocv1_target)
      fullDistrict_vaccinated_ocv2_people <- dplyr::filter(ptargets_avail, fullDistrict_ocv2_target)
      message(paste("Full district ocv1 coverage: ", paste(fullDistrict_vaccinated_ocv1_people$GID_2, collapse = ",")))
      message(paste("Full district ocv2 coverage: ", paste(fullDistrict_vaccinated_ocv2_people$GID_2, collapse = ",")))


      if (nrow(fullDistrict_vaccinated_ocv1_people) == 0){
        message("Partial district coverage: Coverage for ocv1 was so low that no single district was fully vaccinated.")
        partialDistrict_vaccinated_ocv1_people <- goal_ocv1_fvp
        partialDistrict_ocv1_id <- 1

      } else{
        partialDistrict_vaccinated_ocv1_people <- goal_ocv1_fvp - fullDistrict_vaccinated_ocv1_people[which.max(fullDistrict_vaccinated_ocv1_people$possible_fvp_cumsum_ocv1),]$possible_fvp_cumsum_ocv1
        partialDistrict_ocv1_id <- max(fullDistrict_vaccinated_ocv1_people$id)+1
        message(paste("Partial district ocv1 coverage:", ptargets_avail[which(ptargets_avail$id == partialDistrict_ocv1_id),]$GID_2))
      }

      if (nrow(fullDistrict_vaccinated_ocv2_people) == 0){
        message("Partial district coverage: Coverage for ocv2 so low that no single district was fully vaccinated.")
        partialDistrict_vaccinated_ocv2_people <- goal_ocv2_fvp
        partialDistrict_ocv2_id <- 1

      } else{
        partialDistrict_vaccinated_ocv2_people <- goal_ocv2_fvp - fullDistrict_vaccinated_ocv2_people[which.max(fullDistrict_vaccinated_ocv2_people$possible_fvp_cumsum_ocv2),]$possible_fvp_cumsum_ocv2
        partialDistrict_ocv2_id <- max(fullDistrict_vaccinated_ocv2_people$id)+1
        message(paste("Partial district ocv2 coverage:", ptargets_avail[which(ptargets_avail$id == partialDistrict_ocv2_id),]$GID_2))
      }

      ptargets_avail[which(ptargets_avail$id == partialDistrict_ocv1_id),]$actual_ocv1_fvp <- partialDistrict_vaccinated_ocv1_people
      ptargets_avail[which(ptargets_avail$id == partialDistrict_ocv2_id),]$actual_ocv2_fvp <- partialDistrict_vaccinated_ocv2_people


      if(coverage_as_all_0_for_campaign){
        ftargets[[i]] <- dplyr::filter(ptargets_avail, actual_ocv1_fvp>=0 & actual_ocv2_fvp>=0)
      } else if (coverage_as_all_0_for_campaign == FALSE){
        if (scenario == "ocv1-default"){
          ftargets[[i]] <- dplyr::filter(ptargets_avail, actual_ocv1_fvp>0)
        } else {
          ftargets[[i]] <- dplyr::filter(ptargets_avail, actual_ocv1_fvp>0 | actual_ocv2_fvp>0)
        }
      } else{
        stop(message('You did not pass the coverage_as_all_0_for_campaign checkpoint, please go back to the assign_vaccine_targets and check. '))
      }
      names(ftargets)[i] <- as.character(cov_year$year)

    } #endfor

    ftargets_flat <- data.table::rbindlist(ftargets, idcol="vacc_year")

  } #endelse

  ##calam added to export modelled fvps as raw output

  ##set up directory
  #incidence_rate_trend <- as.logical(config$setting$incidence_rate_trend)
  #outbreak_multiplier <- as.logical(config$setting$outbreak_multiplier)
  #setting <- paste0('incid_trend_', incidence_rate_trend, '_outb_layer_',  outbreak_multiplier)
  #dir.create(paste0(rawoutpath, "/","modelled_fvps", "/", scenario, "/", setting), showWarnings = FALSE)
  #mfvps_out_fn <- paste0(rawoutpath, "/","modelled_fvps","/", scenario, "/", setting, "/", country, "_mod_fvps.csv")
  #message(paste("Write modelled fvps:", country, scenario, "\n", ftargets_flat))

  ##write to file
  #readr::write_csv(ftargets_flat, mfvps_out_fn)
  return(ftargets_flat)

}


#' @name run_targeting_strategy
#' @title run_targeting_strategy
#' @description Arrange dataframe so that possible target locations are listed in order of priority
#' @param targets_df dataframe with possible targets
#' @param targeting_strat Unique string to indicate which targeting strategy should be used (currently supports: incidence, affected_pop)
#' @importFrom magrittr %>%
#' @return dataframe that has been rearranged in order that targeting should occur
#' @export
run_targeting_strategy <- function(targets_df, targeting_strat){

  if (targeting_strat == "incidence") {

    if(!"incidence" %in% names(targets_df)){
      stop("Cannot target by incidence when `incidence` is not a column name.")
    }
    rc <- dplyr::arrange(targets_df, desc(incidence))

  } else if (targeting_strat == "random") {   ##order targets randomly for the 'random' targeting strategy
    message("Using random targeting strategy")
    ## vector with random numbers for sorting
    random_id <- sample(1:nrow(targets_df), nrow(targets_df), replace = FALSE)
    targets_df <- targets_df %>%
      dplyr::mutate(sort_id = random_id)
    
    ##order randomly using the random_id column
    rc <- dplyr::arrange(targets_df, sort_id) %>%
      dplyr::select(-sort_id)  ## remove sort_id 
    
  } else if (targeting_strat == "affected_pop"){
    if(!all(c("incidence", "pop_prop") %in% names(targets_df))){
      stop("Cannot target by affected_pop when `incidence` and `pop_prop` are not column names.")
    }
    rc <- dplyr::mutate(targets_df, aff_pop = incidence*pop_prop) %>%
      dplyr::arrange(desc(aff_pop))

  } else{
    stop(paste(targeting_strat, "is not a supported targeting strategy."))
  }

  return(rc)
}




