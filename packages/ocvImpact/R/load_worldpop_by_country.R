
#' @name load_worldpop_by_country
#' @title load_worldpop_by_country
#' @description Load worldpop population raster for a country
#' @param datapath path to data 
#' @param country country code
#' @return population raster 
#' @export
load_worldpop_by_country <- function(datapath, country){

  ## WorldPop population data ##
  pop_fn <- list.files(path = paste0(datapath, "/worldpop"), pattern = "ppp_2020_5km", ignore.case = TRUE)

  if (length(pop_fn)>1){
    stop(paste("There is more than one worldpop file with ppp_2020_5km in the filename"))
  } else if (length(pop_fn)==0){
    stop(paste("There is no worldpop file with  ppp_2020_5km in the filename"))
  } else{
    message(paste0("Loading ", datapath, "/worldpop/", pop_fn))
    pop_world <- raster::raster(paste0(datapath, "/worldpop/", pop_fn))
    pop <- align_rasters(datapath, country, pop_world)
    
    rm(pop_world, shp)
    gc()
  }

  return(pop)
}
