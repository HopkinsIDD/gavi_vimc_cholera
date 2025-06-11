#' @name update_vac_raster
#' @title update_vac_raster
#' @description update_vac_raster one at a time 
#' @param datapath
#' @param modelpath
#' @param country
#' @param scenario
#' @param rawoutpath
#' @param rc_list
#' @param model_year
#' @param pop
#' @param input_list
#' @param no_vacc_year
#' @param rc_targeted
#' @return raster
#' @export
update_vac_raster <- function(datapath,
                              modelpath, 
                              country, 
                              scenario,
                              rawoutpath,
                              rc_list,
                              model_year,
                              pop, # pop raster for the latest year
                              input_list, # input an empty list for the first year, for the following year, input is that list from last year
                              no_vacc_year, 
                              rc_targeted
                              ){

  ## initiate the vac raster 
  message(paste("Stacking population and vaccination proportion raster layer of", model_year, "."))  
  raster0_template <- terra::app(pop, fun = function(x){ifelse(!is.na(x), 0, NA)})
  
  ## quick exit if scenario == 'no-vaccination'
  if(scenario == 'no-vaccination' | no_vacc_year){
    # use the empty raster directly 
    if("rc1" %in% rc_targeted){new_vacc_layer_admin1 <- raster0_template} else {new_vacc_layer_admin1 <- NULL}
    if("rc2" %in% rc_targeted){new_vacc_layer_admin2 <- raster0_template} else {new_vacc_layer_admin2 <- NULL}
    # new_vacc_layer_admin1 <- raster0_template
    # new_vacc_layer_admin2 <- raster0_template

  }else{
    # use the table to guide vaccination 
    if("rc1" %in% rc_targeted){shp1_targeted <- rc_list[1]$rc1[rc_list[1]$rc1$year == model_year & rc_list[1]$rc1$is_target == 1, ]} else {shp1_targeted <- NULL}
    if("rc2" %in% rc_targeted){shp2_targeted <- rc_list[2]$rc2[rc_list[2]$rc2$year == model_year & rc_list[2]$rc2$is_target == 1, ]} else {shp2_targeted <- NULL}
    
    # check point: if there is not a single place targeted this year.
    if(!is.null(shp1_targeted)){

      if(nrow(shp1_targeted) == 0){
        message(paste("There is no admin 1 level areas vaccinated in", model_year, "in", country, "."))
        new_vacc_layer_admin1 <- raster0_template
      }else{
        new_vacc_layer_admin1 <- terra::rasterize(
          terra::vect(shp1_targeted),
          raster0_template,
          field = "actual_prop_vaccinated",
          fun = "mean", # there's one layer which has the same function as "last" in fasterize
          background = 0
        )
        new_vacc_layer_admin1 <- terra::mask(new_vacc_layer_admin1, raster0_template)
      }
      
    }else if(is.null(shp1_targeted)){
      new_vacc_layer_admin1 <- NULL
    }


    if(!is.null(shp2_targeted)){

      if(nrow(shp2_targeted) == 0){
        message(paste("There is no admin 2 level areas vaccinated in", model_year, "in", country, "."))
        new_vacc_layer_admin2 <- raster0_template
      }else{
        new_vacc_layer_admin2 <- terra::rasterize(
          terra::vect(shp2_targeted),
          raster0_template,
          field = "actual_prop_vaccinated",
          fun = "mean", # there's one layer which has the same function as "last" in fasterize
          background = 0
        )
        new_vacc_layer_admin2 <- terra::mask(new_vacc_layer_admin2, raster0_template)
      }
      
    }else if(is.null(shp2_targeted)){
      new_vacc_layer_admin2 <- NULL
    }

  }
  
  ## read in stacks calculated from last campaign year and add new layer or create new ones 
  if(!is.null(input_list)){
    vacc_rasterStack_admin1 <- input_list[["vacc_rasterStack_admin1"]]
    vacc_rasterStack_admin2 <- input_list[["vacc_rasterStack_admin2"]]
    pop_rasterStack <- input_list[["pop_rasterStack"]]
    
    pop_rasterStack <- c(pop_rasterStack, pop)
    if("rc1" %in% rc_targeted){vacc_rasterStack_admin1 <- c(vacc_rasterStack_admin1, new_vacc_layer_admin1)}
    if("rc2" %in% rc_targeted){vacc_rasterStack_admin2 <- c(vacc_rasterStack_admin2, new_vacc_layer_admin2)}
  }else{
    pop_rasterStack <- pop
    vacc_rasterStack_admin1 <- new_vacc_layer_admin1
    vacc_rasterStack_admin2 <- new_vacc_layer_admin2
  }

  ## save to the list 
  input_list = list("vacc_rasterStack_admin1" = vacc_rasterStack_admin1,
                    "vacc_rasterStack_admin2" = vacc_rasterStack_admin2,
                    "pop_rasterStack" = pop_rasterStack)
  rm(new_vacc_layer_admin1, new_vacc_layer_admin2, vacc_rasterStack_admin1, vacc_rasterStack_admin2, pop_rasterStack) 

  return(input_list)
  
}

#' @name save_vac_raster
#' @title save_vac_raster
#' @description save_vac_raster
#' @param datapath
#' @param modelpath
#' @param country
#' @param nsamples
#' @param model_year
#' @param input_list
#' @param rawoutpath
#' @param clean
#' @param scenario
#' @param incidence_rate_trend
#' @param outbreak_multiplier
#' @param vac_incid_threshold
#' @param surveillance_scenario
#' @return save raster
#' @export
save_vac_raster <- function(datapath,
                            modelpath, 
                            country, 
                            nsamples, 
                            model_year,
                            input_list, 
                            rawoutpath = NULL, 
                            clean = NULL, 
                            scenario = NULL, 
                            incidence_rate_trend = NULL, 
                            outbreak_multiplier = NULL, 
                            vac_incid_threshold = NULL, 
                            surveillance_scenario = NULL # input an empty list for the first year, for the following year, input is that list from last year
                            ){
  
  ### Create the file names
  vac1_inter_fn <- paste0("intermediate_raster/", country, "_vac_admin1_", model_year, ".tif")
  vac2_inter_fn <- paste0("intermediate_raster/", country, "_vac_admin2_", model_year, ".tif")
  vac_pop_fn <- paste0("intermediate_raster/", country, "_vac_pop_", model_year, ".tif")

  vac1_out_fn <- paste0(rawoutpath, "/", scenario, "/", paste("incid", incidence_rate_trend, "outbk", outbreak_multiplier, 
                        vac_incid_threshold, surveillance_scenario, country, sep = "_"), "_vac_admin1_", model_year, ".tif")
  vac2_out_fn <- paste0(rawoutpath, "/", scenario, "/", paste("incid", incidence_rate_trend, "outbk", outbreak_multiplier, 
                        vac_incid_threshold, surveillance_scenario, country, sep = "_"), "_vac_admin2_", model_year, ".tif")
  

  ### Convert the structure of the layers 
  if(!is.null(input_list)){
    vac_admin1 <- input_list[[1]]$vacc_rasterStack_admin1
    vac_admin2 <- input_list[[1]]$vacc_rasterStack_admin2
    vac_pop <- input_list[[1]]$pop_rasterStack #one layer is enough for population 

    if(nsamples !=1){
      for(layer_idx in 2:nsamples){
        if(!is.null(vac_admin1)){vac_admin1 <- terra::rast(vac_admin1, input_list[[layer_idx]]$vacc_rasterStack_admin1)}
        if(!is.null(vac_admin2)){vac_admin2 <- terra::rast(vac_admin2, input_list[[layer_idx]]$vacc_rasterStack_admin2)}
      }      
    }

  }else if(!is.null(rawoutpath)){ #this condition applies to when all the vac rasters have been saved in the intermediate folder waiting to be saved in the final output folder
    if( !file.exists(vac1_out_fn) | (file.exists(vac1_out_fn)&clean) ){
      file.rename(from=vac1_inter_fn, to=vac1_out_fn)}
    if( !file.exists(vac2_out_fn) | (file.exists(vac2_out_fn)&clean) ){
      file.rename(from=vac2_inter_fn, to=vac2_out_fn)}
    return(NULL)
  }


  ### Save 
  if(is.null(rawoutpath)){ #if no raw output path is specified, just save the vac rasters in the intermediate folder
    dir.create(paste0("intermediate_raster/"), showWarnings = FALSE)
    if(!is.null(vac_admin1)){terra::writeRaster(vac_admin1, filename = vac1_inter_fn, overwrite = TRUE)}
    if(!is.null(vac_admin2))(terra::writeRaster(vac_admin2, filename = vac2_inter_fn, overwrite = TRUE))
    terra::writeRaster(vac_pop, filename = vac_pop_fn, overwrite = TRUE)
  }else{ #if raw output path is specified, then directly save them in the raw output folder
    message(paste("Writing proportion vaccinated rasterStack for", country))
    dir.create(paste0(rawoutpath, "/", scenario, "/"), showWarnings = FALSE)
    if( !file.exists(vac1_out_fn) | (file.exists(vac1_out_fn)&clean) ){
      if(!is.null(vac_admin1)){terra::writeRaster(vac_admin1, filename = vac1_out_fn, overwrite = TRUE)}}
    if( !file.exists(vac2_out_fn) | (file.exists(vac2_out_fn)&clean) ){
      if(!is.null(vac_admin2)){terra::writeRaster(vac_admin2, filename = vac2_out_fn, overwrite = TRUE)}}
  }

}


