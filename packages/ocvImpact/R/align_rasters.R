#' @name align_rasters
#' @title align_rasters
#' @description Align a raster to the extent and resolution of the worldpop rasters and GADM shapefiles
#' @param datapath path to input data 
#' @param country country code
#' @param orig_raster Raster* object that should be aligned to that in WorldPop and GADM
#' @return Raster* object with aligned extent and resolution at country level
#' @export
#' @include load_shapefile_by_country.R load_worldpop_by_country.R
align_rasters <- function(datapath, country, orig_raster){
  library(terra)
  
  ## if we are using the custom shapefile with health zones (for the DRC case study), specified in the config
  if(as.logical(config$custom$use_custom_shapefile) == TRUE){
    message(paste("Aligning rasters using custom shapefile: ", config$custom$shapefile_filename))
    shp <- load_custom_shapefile_by_country(admin0 = FALSE)
  } else {
    shp <- load_shapefile_by_country(datapath, country, simple=TRUE) ## if we are using the GADM shapefile (VIMC Core model)
    message("Aligning rasters using gadm admin 0 shapefile.")
  }
  pop <- load_worldpop_by_country(datapath, country)
  cropped <- crop(orig_raster, shp, snap = "out")
  masked <- mask(cropped, vect(shp))
  aligned <- resample(masked, rast(pop), method = "near")
  rm(orig_raster, shp, cropped, masked, pop)
  gc()

  return(aligned)
}
