#' @name align_rasters
#' @title align_rasters
#' @description Align a raster to the extent and resolution of the worldpop rasters and GADM shapefiles
#' @param datapath path to input data 
#' @param country country code
#' @param orig_raster Raster* object that should be aligned to that in WorldPop and GADM
#' @return Raster* object with aligned extent and resolution at country level
#' @export
align_rasters <- function(datapath, country, orig_raster){
  
  # shp <- load_shapefile_by_country(datapath, country)
  # masked <- raster::mask(orig_raster, shp, updatevalue = NA)
  # trimmed <- raster::trim(masked)
  shp <- load_shapefile_by_country(datapath, country)
  shp_raster <- raster::rasterize(shp, orig_raster, get_cover = TRUE)
  shp_raster[shp_raster==0] <- NA
  masked <- raster::mask(orig_raster, shp_raster, updatevalue = NA)
  trimmed <- raster::trim(masked)

  rm(orig_raster, masked, shp, shp_raster)
  gc()

  return(trimmed)
}
