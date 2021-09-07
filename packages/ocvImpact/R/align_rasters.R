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
  
  shp <- load_shapefile_by_country(datapath, country, simple=TRUE)
  pop <- load_worldpop_by_country(datapath, country)
  cropped <- raster::crop(orig_raster, shp, snap = "out")
  masked <- raster::mask(cropped, shp, updatevalue = NA) #this is newly added 7/2021
  aligned <- raster::resample(masked, pop, method = "ngb")

  rm(orig_raster, shp, cropped, masked, pop)
  gc()

  return(aligned)
}
