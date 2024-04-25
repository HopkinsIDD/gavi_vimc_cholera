
#' @name load_worldpop_by_country
#' @title load_worldpop_by_country
#' @description Load worldpop population raster for a country
#' @param datapath path to data 
#' @param country country code
#' @return population raster 
#' @export
#' @include load_shapefile_by_country.R
load_worldpop_by_country <- function(datapath, country){

  ## WorldPop population data ##
  pop_fn <- list.files(path = paste0(datapath, "/worldpop"), pattern = "ppp_2020_5km", ignore.case = TRUE)

  if (length(pop_fn)>1){
    stop(paste("There is more than one worldpop file with ppp_2020_5km in the filename"))
  } else if (length(pop_fn)==0){
    stop(paste("There is no worldpop file with  ppp_2020_5km in the filename"))
  } else{
    message(paste0("Loading ", datapath, "/worldpop/", pop_fn))
    
    while(!exists('pop_world')){
      try(pop_world <- raster::raster(paste0(datapath, "/worldpop/", pop_fn)))
      date_time<-Sys.time()
      while((as.numeric(Sys.time()) - as.numeric(date_time))<3.0){}
    }
    # pop_world <- raster::raster(paste0(datapath, "/worldpop/", pop_fn))
    
    ## if we are using the custom shapefile with health zones (DRC case study), specified in config
    if(as.logical(config$use_custom_shapefile) == TRUE){
      message("Use custom shapefile to load worldpop population")
      shp <- load_custom_shapefile_by_country(admin0 = TRUE)
    } else {
      message("Use GADM admin 0 shapefile to load worldpop population")
      shp <- load_shapefile_by_country(datapath, country, simple=TRUE) ## if we are using the GADM shapefile (VIMC Core model)
    }
    cropped <- raster::crop(pop_world, shp, snap = "out")
    pop <- raster::mask(cropped, shp, updatevalue = NA)
    
    rm(pop_world)
    gc()
  }

  return(pop)
}
