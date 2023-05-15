#' @name outbreak_incidence_rate_multiplier
#' @title outbreak_incidence_rate_multiplier
#' @description Create and return a function that can give a prediction of the district-based incidence rate in that country
#' @param datapath path to input data 
#' @param modelpath path to montagu files
#' @param country country code
#' @param lambda the incidence rate raster data that already has been generated
#' @param population_raster the 101-layer population raster that also has already been generated
#' @param output_years the output years for simulation
#' @param use_country_incid_trend whether use the country-level incidence rate trend
#' @param use_original_AR whether directly use the attack rate from the dataset
#' @param draw_district_AR whether just use the attack rate specific to the district
#' @param use_threshold whether use the threshold incidence rate for the non-outbreak months
#' @param smaller_than_1_multiplier_accept whether accept the multipliers that are smaller than 1
#' @param diagnosis whether do diagnosis
#' @param outbreak_df_based_scalar the multiplier/scalar is solely based on the outbreak dataset
#' @return a function that can take in years and give an outbreak multiplier prediction for a given country
#' @export
#' @include utils.R utils_montagu.R
outbreak_incidence_rate_multiplier <- function(datapath,
                                               modelpath,
                                               country,
                                               lambda,
                                               population_raster, 
                                               output_years, 
                                               use_country_incid_trend, 
                                               use_original_AR = TRUE, 
                                               draw_district_AR = TRUE, 
                                               use_threshold = TRUE, 
                                               smaller_than_1_multiplier_accept = FALSE, 
                                               diagnosis = FALSE, 
                                               outbreak_df_based_scalar = TRUE){
  ### Load necessary package
  library(dplyr)
  library(sp)
  library(sf)
  


  ### Get the external outbreak dataset
  outbreak_file_name <- paste0(datapath, '/outbreak/outbreak_df.csv')
  if (file.exists(outbreak_file_name)){
    outbreak_df <- readr::read_csv(outbreak_file_name)
  }else{
    stop('The multi-country outbreak dataset does not exist under the correct directory. Please check. ')
  }
  # What if this dataset does not have data for certain countries 
  if (!country %in% unique(outbreak_df$country_code) | !"admin2" %in% unique(outbreak_df[outbreak_df$country_code == country, ]$spatial_scale)){
    #if there is no outbreak data for the country
    return(function(yr_index){return(1)})
  } else{
    single_country <- outbreak_df %>% filter(country_code == country)
  }
  
  
  
  ### Set the default multiplier raster and make them a list
  multiplier_raster <- lambda
  multiplier_raster[!is.na(multiplier_raster[])] <- 1
  year_length_simulation <- length(output_years)
  multiplier_raster_list <- vector(mode = "list", length = year_length_simulation)
  for (i in 1:year_length_simulation){
    multiplier_raster_list[[i]] <- multiplier_raster
  }
  


  ### Get the year length for the outbreak probability calculation
  year_length_dataset = abs(lubridate::time_length(difftime(min(outbreak_df$TL, na.rm = TRUE), 
                                                            max(outbreak_df$TR, na.rm = TRUE)), "years"))
  day_length_dataset = abs(lubridate::time_length(difftime(min(outbreak_df$TL, na.rm = TRUE), 
                                                           max(outbreak_df$TR, na.rm = TRUE)), "days"))
  


  ### Get the boundary sf
  country2_bd_sf <- rgeoboundaries::gb_adm2(country)


  
  ### Transform the admin2 districts' names
  country2_bd_sf <- country2_bd_sf %>% 
    mutate( NameCode = tolower(gsub("[^a-zA-Z0-9]", "", iconv(shapeName, to = 'ASCII//TRANSLIT'))) )
  single_country_admin2 <- single_country %>% 
    filter(spatial_scale == 'admin2') %>% #we just need the admin2 data
    mutate( Admin2Code = tolower(gsub("[^a-zA-Z0-9]", "", iconv(admin2, to = 'ASCII//TRANSLIT'))) )
  


  ### Get the distribution of AR -- still just based on the admin2-level AR for now (this may need population_raster for AR calculation)
  if(use_original_AR | outbreak_df_based_scalar){
    # Only look into matched locations 
    matched_code_vector <- c()
    for (NameCode in unique(country2_bd_sf$NameCode)){
      if (any( grepl(NameCode, unique(single_country_admin2$Admin2Code), fixed = TRUE) )){
        if (sum( grepl(NameCode, unique(single_country_admin2$Admin2Code), fixed = TRUE) ) > 1){
          ### If more than one match
          matched_code_list <- unique(single_country_admin2[grepl(NameCode, single_country_admin2$Admin2Code, fixed = TRUE), ]$Admin2Code)
          matched_code <- matched_code_list[match(min(nchar(matched_code_list)), nchar(matched_code_list))[1]]
          if (length(match(min(nchar(matched_code_list)), nchar(matched_code_list))) > 1){
            warning('***More than one matched admin2 were identified, the first one was used by default. This warning is from the outbreak_incidence_rate_multiplier function. ')
          }
          matched_code_vector <- c(matched_code_vector, matched_code)
        } else{
          ### If only one match
          matched_code <- unique(single_country_admin2[grepl(NameCode, single_country_admin2$Admin2Code, fixed = TRUE), ]$Admin2Code)
          matched_code_vector <- c(matched_code_vector, matched_code)
        }
      }
    }
    matched_code_vector <- unique(matched_code_vector)

    if(outbreak_df_based_scalar){
      outbreak_scalar_df <- single_country_admin2 %>% 
        filter(Admin2Code %in% matched_code_vector) %>%
        mutate(AR = attack_rate / 1000) %>% #different unit (per 1000 people)
        mutate( middle_date = TL + (abs(difftime(TL, TR, units = "days"))/2) ) %>% 
        mutate( year_new = as.integer(lubridate::year(middle_date)) ) %>% 
        mutate(threshold_new = threshold / 100000) %>% #different unit (per 100,000 people)
        mutate(threshold = threshold_new) %>%
        group_by(Admin2Code) %>% 
        mutate(num_outbreak = length(unique(outbreak_number))) %>% 
        mutate(across_year_AR = (sum(AR * duration) + threshold * (day_length_dataset - sum(duration)))/day_length_dataset) %>% 
        ungroup() %>% 
        rowwise() %>% 
        mutate(outbreak_scalar =  ( (AR * duration + threshold * ((366 - (year_new %% 4)) - duration))/(366 - (year_new %% 4)) ) / across_year_AR) %>% 
        # mutate(outbreak_scalar = AR /  ( (366 - (year_new %% 4))/duration )  / across_year_AR) %>% 
        mutate(threshold_scalar = threshold / across_year_AR) %>% 
        select(Admin2Code, year_new, AR, duration, threshold, num_outbreak, across_year_AR, outbreak_scalar, threshold_scalar)
      
      if(diagnosis){
        write.csv(outbreak_scalar_df, paste0(datapath, '/outbreak/', country, '_outbreak_diag.csv'), row.names = FALSE)
      }
      
    }else{
      AR_duration_df <- single_country_admin2 %>% 
        filter(Admin2Code %in% matched_code_vector) %>% 
        mutate(AR = attack_rate / 1000) %>% #different unit (per 1000 people)
        mutate( middle_date = TL + abs(difftime(TL, TR, units = "days")) ) %>% 
        mutate( year_new = as.integer(lubridate::year(middle_date)) ) %>% 
        mutate(duration_new = duration / (365 + as.integer(year_new%%4 == 0))) %>% 
        mutate(duration = duration_new) %>% 
        mutate(threshold_new = threshold / 100000) %>% #different unit (per 100,000 people)
        mutate(threshold = threshold_new) %>% 
        select(Admin2Code, AR, duration, threshold)
    }
    

  } else {
    AR_vector <- c()
    AR_duration_df <- tibble::tibble(Admin2Code = NA, AR = NA, duration = NA, threshold = NA)
    
    for (NameCode in unique(country2_bd_sf$NameCode)){
      if (any( grepl(NameCode, unique(single_country_admin2$Admin2Code), fixed = TRUE) )){
        if (sum( grepl(NameCode, unique(single_country_admin2$Admin2Code), fixed = TRUE) ) > 1){
          ### If more than one match
          matched_code_list <- unique(single_country_admin2[grepl(NameCode, single_country_admin2$Admin2Code, fixed = TRUE), ]$Admin2Code)
          matched_code <- matched_code_list[match(min(nchar(matched_code_list)), nchar(matched_code_list))[1]]
          if (length(match(min(nchar(matched_code_list)), nchar(matched_code_list))) > 1){
            warning('***More than one matched admin2 were identified, the first one was used by default. This warning is from the outbreak_incidence_rate_multiplier function. ')
          }
        } else{
          ### If only one match
          matched_code <- unique(single_country_admin2[grepl(NameCode, single_country_admin2$Admin2Code, fixed = TRUE), ]$Admin2Code)
        }
        
        ### Get the years of outbreak for a certain district
        single_district_outbreak_ds <- subset(single_country_admin2, Admin2Code == matched_code)
        sCh_vector <- single_district_outbreak_ds %>% 
          mutate( middle_date = TL + abs(difftime(TL, TR, units = "days")) ) %>% 
          mutate( year_new = as.integer(format(middle_date, '%Y')) )
        
        ### Get the population size for that district by year
        single_district_sf <- country2_bd_sf[country2_bd_sf$NameCode == NameCode, ]
        sCh_vector <- sCh_vector %>% arrange(year_new)
        district_pop_vector <- tibble::tibble(year_new = unique(sCh_vector$year_new), 
                                              pop = NA)
        for (real_year in district_pop_vector$year_new){
          yr_index <- real_year - min(output_years) + 1
          dist_population_raster <- raster::subset(population_raster, yr_index, drop = FALSE)
          dist_population_raster <- raster::crop(dist_population_raster, single_district_sf)
          dist_population_raster <- raster::mask(dist_population_raster, single_district_sf)
          pop <- sum(raster::getValues(dist_population_raster), na.rm = TRUE)
          if (pop < 1) {
            warning(paste0('***The total population for the admin2 district ', NameCode, ' is smaller than 1, 1 is used to finish the run. '))
            warning('This above warning is from the outbreak_incidence_rate_multiplier function. ')
            pop <- 1
          }
          district_pop_vector[district_pop_vector$year_new == real_year, ]$pop <- pop
          
        }
        
        ### Join two datasets and output a single simple dataframe
        sCh_pop_vector <- left_join(sCh_vector, district_pop_vector, by = 'year_new') %>% 
          mutate(AR = total_suspected_cases / pop) %>% 
          mutate(duration_new = duration / (365 + as.integer(year_new%%4 == 0))) %>% 
          mutate(duration = duration_new) %>% 
          mutate(threshold_new = threshold / 100000) %>% #different unit (per 100,000 people)
          mutate(threshold = threshold_new) %>% 
          select(Admin2Code, AR, duration, threshold)
        
        AR_duration_df <- rbind(AR_duration_df, sCh_pop_vector)
      }
    }
    #AR_duration_df <- AR_duration_df[-1, ] #get rid of the first row of just NA
    all_na_row <- c(1:nrow(AR_duration_df))[  is.na(AR_duration_df$Admin2Code) 
                                            & is.na(AR_duration_df$AR)
                                            & is.na(AR_duration_df$duration)
                                            & is.na(AR_duration_df$threshold)  ]
    AR_duration_df <- AR_duration_df[ - all_na_row, ] #get rid of the first row of just NA
  }
  
  if(!outbreak_df_based_scalar){
    # Check the emptyness in the AR dataframe
    if( is.null(AR_duration_df$AR) | sum(!is.na(AR_duration_df$AR))==0 ){
      warning('***The AR vector is empty, please check the outbreak dataset for country ', country)
      warning('Single value 1 is returned as the outbreak multiplier for this country. ')
      return(function(yr_index){return(1)})
    } 
  
    # Check NAs in this new simple dataframe
    if (sum(is.na(AR_duration_df)) > 0){
      warning(paste0('***The number of NA(s) in the AR_duration_df is ', sum(is.na(AR_duration_df)), ' in country ', country))
      print(AR_duration_df)
      warning('Now trying to get rid of the NA(s). ')
      AR_duration_df <- na.omit(AR_duration_df)
    }
    
    # Change the name of the new dataframe
    AR_duration_country_df <- AR_duration_df
    rm(AR_duration_df)
  }
  
  

  # #############========= Prepare the incidence rate trend function =========#############
  # incidence_rate_trend <- as.logical(config$setting$incidence_rate_trend)
  # if(incidence_rate_trend){
  #   incid_trend_function <- ocvImpact::generate_flatline_multiplier(
  #                                       trendtype = 'incidence rate', 
  #                                       datapath = datapath, 
  #                                       modelpath = modelpath, 
  #                                       country = country, 
  #                                       use_country_incid_trend = as.logical(use_country_incid_trend))
  # }else{ 
  #   incid_trend_function <- function(year){return(1)}
  # }
  
  #############========= Diagnosis =========#############
  outbreak_diag <- tibble::tibble(NameCode = NA, Admin2Code = NA, outbreak_pro = NA, outbreak_year = NA, 
                                  pop = NA, AR = NA, lambda_incid = NA, outbreak_scalar = NA, 
                                  threshold = NA, non_outbreak_scalar = NA, outbreak_scalar_adjusted = NA)
  
  ### Make an external dataset that saves the indices of the districts
  district_index_fn <- paste0(datapath, "/outbreak/", country, "_district_index.csv")
  if(file.exists(district_index_fn)){
    district_index <- readr::read_csv(district_index_fn)
  }else{
    district_index <- tibble::tibble(DistrictCode = matched_code_vector, row = NA, col = NA)
  }



  ### Go through districts and years
  for (NameCode in unique(country2_bd_sf$NameCode)){
    if (any( grepl(NameCode, unique(single_country_admin2$Admin2Code), fixed = TRUE) )){
      if (sum( grepl(NameCode, unique(single_country_admin2$Admin2Code), fixed = TRUE) ) > 1){
        ### If more than one match
        matched_code_list <- unique(single_country_admin2[grepl(NameCode, single_country_admin2$Admin2Code, fixed = TRUE), ]$Admin2Code)
        matched_code <- matched_code_list[match(min(nchar(matched_code_list)), nchar(matched_code_list))[1]]
        if (length(match(min(nchar(matched_code_list)), nchar(matched_code_list))) > 1){
          warning('***More than one matched admin2 were identified, the first one was used by default. This warning is from the outbreak_incidence_rate_multiplier function. ')
        }
      } else{
        ### If only one match
        matched_code <- unique(single_country_admin2[grepl(NameCode, single_country_admin2$Admin2Code, fixed = TRUE), ]$Admin2Code)
      }

      ### For debugging
      print(paste0("The district being processed now is: ", matched_code))
      
      ### Get the single district sf
      single_district_sf <- country2_bd_sf[country2_bd_sf$NameCode == NameCode, ]
      
      ### What if this district is too small to have any actual values? 
      ### Also record the indices of the current district if not available in the dataset yet
      if(any(is.na( district_index$row ))){
        district_for_NA_diag <- raster::crop(lambda, single_district_sf)
        district_for_NA_diag <- raster::mask(district_for_NA_diag, single_district_sf)

        if(all(is.na(raster::values(district_for_NA_diag))) | is.null(raster::values(district_for_NA_diag))){
          district_index <- district_index %>%
            filter(DistrictCode != matched_code)
          warning(paste0('**The district ', NameCode, ' has been skipped for outbreak multiplier because it is too small to operate. '))
          rm(district_for_NA_diag)
          next
        }else if(any(is.na( district_index[district_index$DistrictCode == matched_code, ]$row ))){
          
          if(dim(district_for_NA_diag)[1] == 1 & dim(district_for_NA_diag)[2] == 1){
            next
          }else{
            district_for_NA_diag <- raster::resample(district_for_NA_diag, lambda) #important
            district_lambda_ext <- raster::extend(district_for_NA_diag, lambda, value = NA) #not sure if necessary
            
            index_matrix <- raster::rowColFromCell(district_lambda_ext, which(!is.na(raster::values(raster::subset(district_lambda_ext, 1)))))
            # index_matrix <- raster::rowColFromCell(lambda, which(!is.na(raster::values(raster::subset(district_lambda_ext, 1)))))
            if( all(is.na(index_matrix)) | is.null(index_matrix) ){
              district_index <- district_index %>%
                filter(DistrictCode != matched_code)
              warning(paste0('**The district ', matched_code, ' does not have indices to assign values, it has been skipped. '))
              next
            }

            district_index_add_on <- district_index[rep( match(matched_code, district_index$DistrictCode), (nrow(index_matrix) - 1) ), ]
            district_index <- rbind(district_index, district_index_add_on) %>% arrange(DistrictCode)
            district_index[district_index$DistrictCode == matched_code, ]$row <- index_matrix[, 1]
            district_index[district_index$DistrictCode == matched_code, ]$col <- index_matrix[, 2]
            rm(district_for_NA_diag)
            rm(district_lambda_ext)
          }
          
        }
      }
      
      
      
      ### Calculate the annual outbreak probability in this district
      single_district_outbreak_ds <- subset(outbreak_scalar_df, Admin2Code == matched_code)
      outbreak_pro <- single_district_outbreak_ds$num_outbreak[1] / year_length_dataset
      
      
      
      ### One layer in population_raster represents one year 
      for (yr_index in 1:year_length_simulation) {
        ## Get the layers that will have outbreak
        outbreak_layer_index <- sample(1:as.numeric(raster::nlayers(lambda)), 
                                      size = as.numeric(raster::nlayers(lambda))*outbreak_pro, 
                                      replace = FALSE)
        if(length(outbreak_layer_index) == 0){
          warning(paste0("*** The district ", matched_code, " will not have an outbreak multiplier because the probability of outbreak for it is too low. "))
          next
        }

        ## Get the outbreak_scalar
        if(length(single_district_outbreak_ds$outbreak_scalar) == 1){
          outbreak_scalar <- rep(single_district_outbreak_ds$outbreak_scalar, length(outbreak_layer_index))
        }else{
          outbreak_scalar <- sample(single_district_outbreak_ds$outbreak_scalar, 
                                    size = length(outbreak_layer_index), replace = TRUE)
        }

        ## Assign outbreak scalars
        multiplier_raster_for_change <- multiplier_raster_list[[yr_index]]
        for(i in outbreak_layer_index){
          for(h in unique(district_index[district_index$DistrictCode == matched_code, ]$row)){
            multiplier_raster_for_change[[i]][h, 
                                              district_index[district_index$DistrictCode == matched_code & 
                                                             district_index$row == h, ]$col ] <- outbreak_scalar[match(i, outbreak_layer_index)]
          }
        }
        
        ## Get the layers that will not have outbreak
        threshold_layer_index <- (1:as.numeric(raster::nlayers(lambda)))[!(1:as.numeric(raster::nlayers(lambda))) %in% outbreak_layer_index]

        ## Get the threshold_layer
        if(length(single_district_outbreak_ds$threshold_scalar) == 1){
          threshold_scalar <- rep(single_district_outbreak_ds$threshold_scalar, length(threshold_layer_index))
        }else{
          threshold_scalar <- sample(single_district_outbreak_ds$threshold_scalar, 
                                    size = length(threshold_layer_index), replace = TRUE)
        }

        ## Assign threshold scalars
        for(j in threshold_layer_index){
          for(k in unique(district_index[district_index$DistrictCode == matched_code, ]$row)){
            multiplier_raster_for_change[[j]][k, 
                                              district_index[district_index$DistrictCode == matched_code & 
                                                             district_index$row == k, ]$col ] <- threshold_scalar[match(j, threshold_layer_index)]
          }
        }
        
        ## Put back the raster
        multiplier_raster_list[[yr_index]] <- multiplier_raster_for_change

      } 
    } 
  } #this is the other side of the for loop 
  


  ### Save the diagnosis dataframe
  if(diagnosis & !outbreak_df_based_scalar){
    outbreak_diag <- outbreak_diag[-1, ] #delete the NA row
    
    outbreak_diag = outbreak_diag %>% 
      tidyr::gather(key = 'VariableName', value = 'VariableValue', AR:outbreak_scalar_adjusted) %>% 
      tidyr::separate(VariableValue, paste0('layer_', as.character(1:as.numeric(raster::nlayers(lambda)))), sep = ",") %>% 
      arrange(NameCode, outbreak_year)
    
    invisible(sapply(outbreak_diag[c('outbreak_pro', 'pop')], as.numeric))
    #invisible(sapply(unlist(outbreak_diag[paste0('layer_', as.character(1:as.numeric(raster::nlayers(lambda))))]), as.numeric))
    for ( var_idx in 1:as.numeric(raster::nlayers(lambda)) ) {
      outbreak_diag[, 6 + var_idx] <- as.numeric(unlist(outbreak_diag[, 6 + var_idx]))
    }
    
    outbreak_diag <- outbreak_diag %>%
      mutate_if(is.numeric, round, digits = 7)
    
    write.csv(outbreak_diag, paste0(datapath, '/outbreak/', country, '_outbreak_diag.csv'), row.names = FALSE)
  }
  
  ### Save the raster file if it only represents one year -- countrycode_settingnumber_year.tif
  if(length(multiplier_raster_list) == 1){
    random_seed <- as.numeric(config$setting$random_seed)
    setting_num <- random_seed #for now, yes
    outbreak_out_fn <- paste0(datapath, '/outbreak/', country, '_', (output_years[1] %% 10), '.tif')
    raster::writeRaster(raster::stack(multiplier_raster_list[[1]]), filename = outbreak_out_fn)
  }

  ### Return the function with a list of multi-layer rasters and save the indices
  if( !file.exists(district_index_fn) ){
    readr::write_csv(district_index, district_index_fn)
  }else if( any(is.na(readr::read_csv(district_index_fn))) ){
    readr::write_csv(district_index, district_index_fn)
  }

  return(function(yr_index){return( multiplier_raster_list[[yr_index]] )})
  
}

# ### tmp: diag
# pdf("/home/kaiyuezou/VIMC_Model/202110gavi-3/gavi_vimc_cholera/test_ben13.pdf", width=10, height=15)
# raster::plot(ben_raster[[1]], useRaster=FALSE)
# raster::plot(ben_raster[[2]], useRaster=FALSE)
# raster::plot(ben_raster[[3]], useRaster=FALSE)
# raster::plot(ben_raster[[4]], useRaster=FALSE)
# raster::plot(ben_raster[[5]], useRaster=FALSE)
# raster::plot(ben_raster[[6]], useRaster=FALSE)
# raster::plot(ben_raster[[7]], useRaster=FALSE)
# raster::plot(ben_raster[[8]], useRaster=FALSE)
# raster::plot(ben_raster[[9]], useRaster=FALSE)
# raster::plot(ben_raster[[10]], useRaster=FALSE)
# dev.off()

