#' @name case_by_district
#' @title case_by_district
#' @description sum number of cases by each district
#' @param case_raster_multilayer the multi-layer expected case raster  
#' @param admin_shp the shapefiles for a certain admin level 
#' @return 
#' @export
case_by_district <- function(case_raster_multilayer, admin_shp){

    true_case_all <- exactextractr::exact_extract(case_raster_multilayer, admin_shp, 'sum')
    return(true_case_all)
}
