
#' @name create_incid_raster
#' @title create_incid_raster
#' @description Generate a raster with the baseline cholera incidence
#' @param datapath path to data 
#' @param country country code
#' @param nsamples numeric, number of layers to sample (must be below 1000)
#' @param clean logical that indicates whether existing vacc_files should be deleted
#' @return raster of incidence rate, 30 samples
#' @export
create_incid_raster <- function(datapath, country, nsamples, clean){

  incid_out_fn <- paste0(datapath, "/", country, "_incid_5k_", nsamples, ".tif")
  
  if(clean & file.exists(incid_out_fn)){
    message(paste("Clean existing", incid_out_fn))
    file.remove(incid_out_fn)
  }

  if(!file.exists(incid_out_fn)){
    if (country %in% c("COD", "ETH", "KEN", "SOM", "SSD")){

      layer_indexes <- sort(sample(1:1000, nsamples, replace=TRUE))
      print(layer_indexes)

      ## incidence data ##
      message(paste0("Loading ", datapath, "/incidence/afro_2010-2016_lambda_5k.tif"))
      afr <- raster::stack(paste0(datapath, "/incidence/afro_2010-2016_lambda_5k.tif"))
      afr_sample <- raster::subset(afr, layer_indexes, drop = FALSE)
      rm(afr)
      gc()

      lambda <- align_rasters(datapath, country, afr_sample)
      rm(afr_sample)
      gc()

      message(paste("Write", incid_out_fn))
      raster::writeRaster(raster::stack(lambda), filename = incid_out_fn)

    } else {
      stop(paste(country, "cholera incidence data was not found."))
    }

  } else {
    message(paste("Skip creation", incid_out_fn))
    lambda <- raster::stack(incid_out_fn)
  }
  
  return(lambda)
}
