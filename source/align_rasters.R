#' @name align_rasters
#' @title align_rasters
#' @description Align a raster to the extent and resolution of the worldpop rasters and GADM shapefiles
#' @param datapath path to input data 
#' @param country country code
#' @param orig_raster Raster* object that should be aligned to that in WorldPop and GADM
#' @return Raster* object with aligned extent and resolution
#' @export
align_rasters <- function(datapath, country, orig_raster, tol = 1e-6){
  pop <- load_worldpop_by_country(datapath, country)
  shp <- load_shapefile_by_country(datapath, country)
  
  masked <- raster::mask(orig_raster, shp, updatevalue = NA)
  trimmed <- raster::trim(masked)
  message("Resampling original raster to match extent of population raster")
  # extended <- raster::extend(trimmed, pop)
  resampled <- raster::resample(extended, pop, method = "bilinear")  

  # if (all(raster::extent(trimmed) != raster::extent(pop))){
  #   message("Resampling original raster to match extent of population raster")
  #   resampled <- raster::resample(trimmed, pop, method = "bilinear")  
  # } else{
  #   message("Original raster matches extent of population raster")
  #   resampled <- trimmed
  # }

  # rc <- fix_raster_tolerance(resampled)

  rm(orig_raster, pop, masked, trimmed, extended)
  gc()

  return(resampled)
}


# #' @name fix_raster_tolerance
# #' @title fix_raster_tolerance
# #' @description Round extent and resolution to meet tolerance
# #' @param orig_raster Raster* object that should be aligned to that in WorldPop and GADM
# #' @return Raster* object with tolerance rounded extent and resolution
# #' @export
# fix_raster_tolerance <- function(orig_raster, tol = 1e-6){
#     raster::extent(orig_raster) <- round(raster::extent(orig_raster), -log10(tol))
#     raster::res(orig_raster) <- round(raster::res(orig_raster), -log10(tol))
#     return(orig_raster)
# }
