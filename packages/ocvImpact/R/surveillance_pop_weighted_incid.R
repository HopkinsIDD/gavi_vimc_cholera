#' @name pop_weighted_admin_mean_incid
#' @title pop_weighted_admin_mean_incid
#' @description pop_weighted_admin_mean_incid
#' @param datapath path to input data 
#' @param modelpath
#' @param incidence_rate_raster
#' @param pop_raster
#' @param country
#' @param admin_shp
#' @return 
#' @export
#' @include
pop_weighted_admin_mean_incid <- function(datapath, modelpath, incidence_rate_raster, pop_raster, country, admin_shp){
  
  pop_for_weights <- raster::resample(pop_raster, incidence_rate_raster)
  case_raster <- incidence_rate_raster * pop_for_weights
  
  # pop_for_weights <- pop_raster
  # case_raster <- raster::overlay(incidence_rate_raster, pop_for_weights, fun = function(x, y){ x*y },
  #                               recycle = TRUE, unstack = TRUE) 
  # case_vector <- raster::extract(case_raster, admin_shp, method = 'bilinear', fun = sum, na.rm = TRUE)
  # pop_vector <- raster::extract(pop_for_weights, admin_shp, method = 'bilinear', fun = sum, na.rm = TRUE)
  case_vector <- exactextractr::exact_extract(case_raster, admin_shp, 'sum')
  pop_vector <- exactextractr::exact_extract(pop_for_weights, admin_shp, 'sum')
  incid_vector <- case_vector / pop_vector

  # case_raster <- stars::st_as_stars(case_raster)
  # case_vector <- nngeo::raster_extract(case_raster, admin_shp, fun = sum, na.rm = TRUE)
  # incid <- exactextractr::exact_extract(incidence_rate_raster, admin_shp, 'weighted_mean', weights = area(pop_for_weights))
  
  return(incid_vector)
}
