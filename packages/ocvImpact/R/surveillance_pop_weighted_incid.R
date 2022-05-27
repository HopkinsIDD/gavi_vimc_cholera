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
  incid <- exactextractr::exact_extract(incidence_rate_raster, admin_shp, 'weighted_mean', weights = area(pop_for_weights))
  
  return(incid)
}
