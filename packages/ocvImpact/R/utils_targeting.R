
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
#' @param country country code
#' @importFrom magrittr %>%
#' @return dataframe with incidence and population by admin unit for a single country
#' @export 
#' @include load_worldpop_by_country.R load_shapefile_by_country.R utils_targeting.R
load_targets_by_country <- function(datapath, country){

  if (country %in% c("COD", "ETH", "KEN", "SOM", "SSD")){

    ## incidence data ##
    message(paste0("Loading ", datapath, "/incidence/afro_2010-2016_lambda_5k_mean.tif"))
    afr <- raster::raster(paste0(datapath, "/incidence/afro_2010-2016_lambda_5k_mean.tif"))

    ## WorldPop population data ##
    pop <- load_worldpop_by_country(datapath, country)
    ## admin unit shapefile ##
    shp <- load_shapefile_by_country(datapath, country)

    ## summarize rasters to admin level
    incid2 <- exactextractr::exact_extract(afr, shp, 'mean') ## could add population weight here for better incidence estimate but need to project population to the incidence grid
    pop2 <- get_admin_population(pop, shp)
    total_pop <- sum(pop2)

    rc <- dplyr::mutate(shp, 
                        incidence = incid2,
                        pop_wp = pop2,
                        pop_prop = pop2/total_pop) %>%
      sf::st_drop_geometry() %>%
      dplyr::select(GID_0, GID_2, NAME_1, NAME_2, incidence, pop_prop) %>%
      tibble::as_tibble()

    if (sum(rc$pop_prop)!=1){
      stop(paste("The population proportion calculation is incorrect for", country))
    }

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
#' @param targeting_strat Unique string that identifies how targets should be ranked (default = incidence)
#' @param campaign_cov Proportion of the target population vaccinated by the campaign (default = 0.8)
#' @param num_skip_years Skip targeting in locations vaccinated within the last X years (default = 0 --> no skipping)
#' @importFrom magrittr %>%
#' @return 
#' @export
#' @include utils_targeting.R utils_montagu.R
assign_vaccine_targets <- function(datapath, modelpath, country, scenario, targeting_strat = "incidence", campaign_cov = 0.8, num_skip_years = 0){

  message(paste("Now assigning vaccine by incidence:", country, scenario))
  ptargets <- load_targets_by_country(datapath, country)
  coverage <- import_coverage_scenario(modelpath, country, scenario, filter0 = TRUE) 

  if (is.null(coverage)){

    message("No vaccine was assigned for the no-vaccination scenario. Returning NULL value.")
    ftargets_flat <- NULL
  
  } else{

    ## Perform checks on the coverage scenario
    if (!all(coverage$gender == "both") | 
        !all(coverage$age_range_verbatim == "default age and gender") | 
        !all(coverage$activity_type == "campaign")
        ){
      stop(paste("Vaccine assignment is not supported for this coverage scenario. Check the gender, age_range_verbatim, and activity_type columns in the", scenario, "coverage sheet."))
    }

    ftargets <- vector(mode = "list", length = length(coverage$year))

    for (i in 1:length(coverage$year)){

      cov_year <- coverage[i,]
      goal_target_pop <- cov_year$target
      goal_fvp <- cov_year$fvp
      
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
          possible_fvp = round(goal_target_pop * pop_prop * campaign_cov, 0), 
          possible_fvp_cumsum = cumsum(possible_fvp),
          actual_fvp = ifelse(possible_fvp_cumsum <= goal_fvp, possible_fvp, 0),
          fullDistrict_target = ifelse(actual_fvp > 0, TRUE, FALSE))

      ## partial district targeting
      fullDistrict_vaccinated_people <- dplyr::filter(ptargets_avail, fullDistrict_target) 
      message(paste("Full district coverage: ", paste(fullDistrict_vaccinated_people$GID_2, collapse = ",")))

      if (nrow(fullDistrict_vaccinated_people) == 0){
        message("Partial district coverage: Coverage was  so low that no single district was fully vaccinated.")
        partialDistrict_vaccinated_people <- goal_fvp
        partialDistrict_id <- 1
      
      } else{
        partialDistrict_vaccinated_people <- goal_fvp - fullDistrict_vaccinated_people[which.max(fullDistrict_vaccinated_people$possible_fvp_cumsum),]$possible_fvp_cumsum
        partialDistrict_id <- max(fullDistrict_vaccinated_people$id)+1
        message(paste("Partial district coverage:", ptargets_avail[which(ptargets_avail$id == partialDistrict_id),]$GID_2))
      }
      
      ptargets_avail[which(ptargets_avail$id == partialDistrict_id),]$actual_fvp <- partialDistrict_vaccinated_people

      ftargets[[i]] <- dplyr::filter(ptargets_avail, actual_fvp>0)
      names(ftargets)[i] <- as.character(cov_year$year)

    } #endfor

    ftargets_flat <- data.table::rbindlist(ftargets, idcol="vacc_year")

  } #endelse

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

      


