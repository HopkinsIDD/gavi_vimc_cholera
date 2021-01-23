#' @name align_rasters
#' @title align_rasters
#' @description Align a raster to the extent and resolution of the worldpop rasters and GADM shapefiles
#' @param datapath path to input data 
#' @param country country code
#' @param orig_raster Raster* object that should be aligned to that in WorldPop and GADM
#' @return Raster* object with aligned extent and resolution at country level
#' @export
align_rasters <- function(datapath, country, orig_raster){
  
  shp <- load_shapefile_by_country(datapath, country, simple=TRUE)
  cropped <- raster::crop(orig_raster, shp, snap = "out")

  rm(orig_raster, shp)
  gc()

  return(cropped)
}
