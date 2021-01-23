
#' @name load_shapefile_by_country
#' @title load_shapefile_by_country
#' @description Load admin2 level shapefile for a country
#' @param datapath path to data 
#' @param country country code
#' @param simple logical indicating whether simple (full country) shapefile should be used (default: FALSE)
#' @return shapefile in sf format
#' @export 
load_shapefile_by_country <- function(datapath, country, simple = FALSE){

  if (country %in% c("COD", "ETH", "KEN", "SOM", "SSD")){
    if (simple){
      country_pattern <- paste(country, "0", sep = "_")
    } else{
      country_pattern <- paste(country, "2", sep = "_")
    }

    shp_fn <- list.files(path = paste0(datapath, "/shapefiles"), pattern = country_pattern)
    if (length(shp_fn) > 1){
      stop(paste("There is more than one shapefile with", country_pattern, "in the filename"))
    }
    message(paste0("Loading ", datapath, "/shapefiles/", shp_fn))
    shp <- readRDS(paste0(datapath, "/shapefiles/", shp_fn))
  
  } else{
    stop(paste(country, "is not yet supported."))
  }

  return(shp)
}
