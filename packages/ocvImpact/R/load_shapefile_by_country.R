
#' @name load_shapefile_by_country
#' @title load_shapefile_by_country
#' @description Download admin2 or admin0 (simple) level shapefile for a country from GADM, if not already in shapefiles folder, and return sf object
#' @param datapath path to data 
#' @param country country code
#' @param simple logical indicating whether simple (full country) shapefile should be used (default: FALSE)
#' @return shapefile in sf format
#' @export 
load_shapefile_by_country <- function(datapath, country, simple = FALSE){

  
  tryCatch(
    {
      if (simple){
          country_pattern <- paste(country, "0", sep = "_")
          #shp <- GADMTools::gadm_sf_loadCountries(c(country), level = 0, basefile = file.path(datapath, "shapefiles/"))$sf
          shp <- geodata::gadm(c(country), level = 0, path = file.path(datapath, "shapefiles/")) ## new implementation 1 Aug 2024
      } else{
          country_pattern <- paste(country, "2", sep = "_")
          #shp <- GADMTools::gadm_sf_loadCountries(c(country), level = 2, basefile = file.path(datapath, "shapefiles/"))$sf
          shp <- geodata::gadm(c(country), level = 2, path =  file.path(datapath, "shapefiles/")) ## new implementation 1 Aug 2024
      }
      ## new implementation 1 Aug 2024
      shp <- sf::st_as_sf(shp, crs = sf::st_crs(4326)) 
      shp <- sf::st_cast(shp, "POLYGON") ## to ensure all rows are of the same geometry and avoid fasterize errors
      message(paste0("Loading ", datapath, "/shapefiles/", country_pattern, "_sf.rds"))
    },
    error = function(e) {
      print(paste("Unable to get shapefile for", country, ".", e))
    }
  )
  ## make shapefile valid, this may fail and throw error on idmodeling because lwgeom is not installed there
  shp <- sf::st_make_valid(shp)
  return(shp)
}
