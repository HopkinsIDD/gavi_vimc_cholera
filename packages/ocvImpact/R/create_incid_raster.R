
#' @name create_incid_raster
#' @title create_incid_raster
#' @description Generate a raster with the baseline cholera incidence
#' @param modelpath path to montagu 
#' @param datapath path to data 
#' @param country country code
#' @param nsamples numeric, number of layers to sample (must be below 1000)
#' @param redraw logical that indicates whether existing incid samples should be drawn again
#' @param random_seed random seed number
#' @return raster of incidence rate, "nsamples" samples
#' @export
#' @include get_singular_estimate.R align_rasters.R utils_montagu.R load_worldpop_by_country.R
#######Kaiyue Added on 7/12/2021#######
###We also include utils_montagu.R just to get population estimate for a country
###We also include load_worldpop_by_country.R just to get a proper raster
###We also include modelpath as input
###########Comment completed###########

create_incid_raster <- function(modelpath, datapath, country, nsamples, redraw, random_seed = NULL){

  #######Kaiyue Added on 7/15/2021#######
  ######Kaiyue editted on 1/30/2022######
  ###Use the WHO data source to update the current spreadsheet
  IncidenceTable <- ocvImpact::get_singular_estimate(datapath, country)
  RasterCountry <- subset(IncidenceTable, is.na(IncidenceTable$singular_estimate))$country_code
  NonRasterCountry <- subset(IncidenceTable, !is.na(IncidenceTable$singular_estimate))$country_code
  ###Calculate and generate
  if (country %in% NonRasterCountry){
    NRCountryIndex <- match(country, IncidenceTable$country_code)
    YearList <- as.numeric(strsplit(IncidenceTable$year_list[NRCountryIndex], '-')[[1]])
    CountryPopTable <- ocvImpact::import_country_population(modelpath, country, redownload = FALSE)
    CountryPopMean <- mean(CountryPopTable$pop_model[match(YearList, CountryPopTable$year)], na.rm = TRUE)
    
    if (IncidenceTable$is_num_case[NRCountryIndex] == 0){
      NumCases <- IncidenceTable$singular_estimate[NRCountryIndex] * CountryPopMean
      IncidenceTable$num_case[NRCountryIndex] <- NumCases
      IncidenceTable$incid_rate_100k[NRCountryIndex] <- 100000*NumCases/CountryPopMean
    } else if (IncidenceTable$is_num_case[NRCountryIndex] == 1){
      NumCases <- IncidenceTable$singular_estimate[NRCountryIndex]
      IncidenceTable$num_case[NRCountryIndex] <- NumCases
      IncidenceTable$incid_rate_100k[NRCountryIndex] <- 100000*NumCases/CountryPopMean
    } else{
      message('The "is_num_case" variable has an invalid value. ')
    }
    
    if(country == "IND"){ #borrow BGD incidence rate for IND
      NumCases <- 0.0008536926*CountryPopMean
    }

    PoisCases <- rpois(nsamples, NumCases)
    StochasticIR <- PoisCases / CountryPopMean
    write.csv(IncidenceTable, file = paste0(datapath, '/incidence/VIMC-47-countries-for-cholera-modelling.csv'), row.names=FALSE)
    rm(IncidenceTable)
    rm(CountryPopTable)
    rm(PoisCases)
    gc()
  }
  ###########Comment completed###########
  
  incid_out_fn <- paste0(datapath, "/incidence/", country, "_incid_5k_", nsamples, ".tif")
  
  if (!file.exists(incid_out_fn) | redraw){
    
    ###If we want to redraw, first delete
    if(file.exists(incid_out_fn)){
      message(paste("Clean existing", incid_out_fn))
      file.remove(incid_out_fn)
    }

    ###For BGD
    if(country == "BGD"){
      if(is.null(random_seed)){
        random_seed <- as.numeric(config$setting$random_seed)
      }
      setting_num <- random_seed #for the current design
      BGD_raster <- raster::stack(paste0(datapath, "/incidence/BGD_incid_5k_50_setting", setting_num, ".tif"))
      return(BGD_raster)
    }
    


    ###Just generate
    #######Kaiyue Added on 7/12/2021#######
    #if (country %in% c("COD", "ETH", "KEN", "SOM", "SSD")){ ----the older code with just the testing countries
    #the next line is the new code
    if (country %in% RasterCountry){
      ###########Comment completed###########
      if(is.null(random_seed)){
        random_seed <- as.numeric(config$setting$random_seed)
      }
      set.seed(random_seed)
      layer_indexes <- sort(sample(1:1000, nsamples, replace=TRUE))
      print(layer_indexes)
      
      ## incidence data ##
      ##addition to use new mai rate raster for 2016-2020 for new touchstone
      runname <- config$runname
      if (runname == "202310gavi-4"){
        message(paste0("Loading ", datapath, "/incidence/afro_2016-2020_lambda_5k.tif"))
        afr <- raster::stack(paste0(datapath, "/incidence/afro_2016-2020_lambda_5k.tif"))
      } else {
        message(paste0("Loading ", datapath, "/incidence/afro_2010-2016_lambda_5k.tif"))
        afr <- raster::stack(paste0(datapath, "/incidence/afro_2010-2016_lambda_5k.tif"))
      }
      ##end addition
      ###########this is following line is just for temp use before we figure out how to transfer raster files without corruptions###########
      #afr <- raster::stack(paste0("/home/kaiyuezou/montagu_try_7_5/gavi_vimc_cholera/input_data/incidence/afro_2010-2016_lambda_5k.tif"))
      afr_sample <- raster::subset(afr, layer_indexes, drop = FALSE)
      rm(afr)
      gc()
      
      lambda <- align_rasters(datapath, country, afr_sample)
      rm(afr_sample)
      gc()
      
      message(paste("Write", incid_out_fn))
      raster::writeRaster(raster::stack(lambda), filename = incid_out_fn)
      
      #######Kaiyue Added on 7/12/2021####### -- this is to use singular incidence to represent the whole country incidence
   }else if(country %in% NonRasterCountry){
      ###Calculate the number of cases for each grid cell and stack
      CountryPopRaster <- ocvImpact::load_worldpop_by_country(datapath, country)
      StochasticLayers <- list()
      for (i in 1:length(StochasticIR)){
        CountryIRRaster <- CountryPopRaster
        CountryIRRaster[!is.na(CountryIRRaster)] <- StochasticIR[i]
        StochasticLayers[[i]] <- CountryIRRaster
      }
      rm(CountryPopRaster)
      rm(CountryIRRaster)
      
      nrc <- raster::stack(StochasticLayers)
      nrc_sample <- nrc #we already have "nsample" samples 
      #raster::plot(nrc_sample) ###################just for testing
      
      rm(StochasticLayers)
      rm(nrc)
      gc()
      lambda <- align_rasters(datapath, country, nrc_sample)
      #raster::plot(lambda) ###################just for testing
      rm(nrc_sample)
      gc()
      message(paste("Write", incid_out_fn))
      raster::writeRaster(raster::stack(lambda), filename = incid_out_fn)
      ###########Comment completed###########
      
   }else {
      stop(paste(country, "cholera incidence data was not found."))
    }
    
 }else {
    message(paste("Skip creation", incid_out_fn))
    lambda <- raster::stack(incid_out_fn)
 }
  
  return(lambda)
}
