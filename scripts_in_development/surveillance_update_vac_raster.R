# Function: update_input_rasterStack() to get the input rasterstack
update_vac_raster <- function(datapath,
                              modelpath, 
                              country, 
                              scenario,
                              rawoutpath,
                              rc_list,
                              model_year,
                              pop, # pop raster for the latest year
                              input_list # input an empty list for the first year, for the following year, input is that list from last year
                              ){

  # initiate the vac raster 
  message(paste("Stacking population and vaccination proportion raster layer of", model_year, "."))  
  raster0_template <- raster::calc(pop, fun = function(x){ifelse(!is.na(x), 0, NA)})
  
    
  # use the table to guide vaccination 
  shp1_targeted <- rc_list[[1]][rc_list[[1]]$year == model_year & rc_list[[1]]$is_target == 1, ]
  shp2_targeted <- rc_list[[2]][rc_list[[2]]$year == model_year & rc_list[[2]]$is_target == 1, ]
  
  # check point: if there is not a single place targeted this year.
  if(nrow(shp1_targeted) == 0){
    
    message(paste("There is no admin 1 level areas vaccinated in", model_year, "in", country, "."))
    new_vacc_layer_admin1 <- raster0_template
  
  }else{
    new_vacc_layer_admin1 <- fasterize::fasterize(
      shp1_targeted,
      raster0_template,
      field = "actual_prop_vaccinated",
      fun = "last",
      background = 0
    )
    new_vacc_layer_admin1 <- raster::mask(new_vacc_layer_admin1, raster0_template, updatevalue = NA)
  }

  if(nrow(shp2_targeted) == 0){
    
    message(paste("There is no admin 2 level areas vaccinated in", model_year, "in", country, "."))
    new_vacc_layer_admin2 <- raster0_template
  
  } else{
    new_vacc_layer_admin2 <- fasterize::fasterize(
      shp2_targeted,
      raster0_template,
      field = "actual_prop_vaccinated",
      fun = "last",
      background = 0
    )
    new_vacc_layer_admin2 <- raster::mask(new_vacc_layer_admin2, raster0_template, updatevalue = NA)
  }
  
  # read in stacks calculated from last campaign year and add new layer or create new ones 
  if(!is.null(input_list)){
    vacc_rasterStack_admin1 <- input_list[["vacc_rasterStack_admin1"]]
    vacc_rasterStack_admin2 <- input_list[["vacc_rasterStack_admin2"]]
    pop_rasterStack <- input_list[["pop_rasterStack"]]
    
    pop_rasterStack <- raster::addLayer(pop_rasterStack, pop)
    vacc_rasterStack_admin1 <- raster::addLayer(vacc_rasterStack_admin1, new_vacc_layer_admin1)
    vacc_rasterStack_admin2 <- raster::addLayer(vacc_rasterStack_admin2, new_vacc_layer_admin2)
  }else{
    pop_rasterStack <- pop
    vacc_rasterStack_admin1 <- new_vacc_layer_admin1
    vacc_rasterStack_admin2 <- new_vacc_layer_admin2
  }

  # save to the list 
  input_list = list("vacc_rasterStack_admin1" = vacc_rasterStack_admin1,
                    "vacc_rasterStack_admin2" = vacc_rasterStack_admin2,
                    "pop_rasterStack" = pop_rasterStack)
  rm(new_vacc_layer_admin1, new_vacc_layer_admin2, vacc_rasterStack_admin1, vacc_rasterStack_admin2, pop_rasterStack) 

  return(input_list)
  
}



### This function will not be called because the first one will suffice from now 
create_first_year_vac_raster <- function( datapath, modelpath, country,
                                          model_year, rc_list,
                                          pop){

    message("Loading first layer of population and vaccination proportion rasterStack.")
    
    raster0_template <- raster::calc(pop, fun = function(x){ifelse(!is.na(x), 0, NA)})
    shp1_targeted <- rc_list[[1]][rc_list[[1]]$year == model_year & rc_list[[1]]$is_target == 1, ]
    shp2_targeted <- rc_list[[2]][rc_list[[2]]$year == model_year & rc_list[[2]]$is_target == 1, ]

    # check point: if there is not a single place targeted this year.
    if(nrow(shp1_targeted) == 0){
      
      message(paste("There is no admin 1 level areas vaccinated in", model_year, "in", country, "."))
      startvacc_raster_admin1 <- raster0_template
    
    }else{
      startvacc_raster_admin1 <- fasterize::fasterize(
        shp1_targeted,
        raster0_template,
        field = "actual_prop_vaccinated",
        fun = "last",
        background = 0
      )
      startvacc_raster_admin1 <- raster::mask(startvacc_raster_admin1, raster0_template, updatevalue = NA)
    }

    if(nrow(shp2_targeted) == 0){
      
      message(paste("There is no admin 2 level areas vaccinated in", model_year, "in", country, "."))
      startvacc_raster_admin2 <- raster0_template
    
    }else{
      startvacc_raster_admin2 <- fasterize::fasterize(
        shp2_targeted,
        raster0_template,
        field = "actual_prop_vaccinated",
        fun = "last",
        background = 0
      )
      startvacc_raster_admin2 <- raster::mask(startvacc_raster_admin2, raster0_template, updatevalue = NA)
    }

    # stack
    vacc_rasterStack_admin1 <- raster::stack(startvacc_raster_admin1)
    vacc_rasterStack_admin2 <- raster::stack(startvacc_raster_admin2)
    pop_rasterStack <- raster::stack(pop)

    input_list = list("vacc_rasterStack_admin1" = vacc_rasterStack_admin1,
                      "vacc_rasterStack_admin2" = vacc_rasterStack_admin2,
                      "pop_rasterStack" = pop_rasterStack)
    
    rm(vacc_rasterStack_admin1, vacc_rasterStack_admin2, pop_rasterStack)

    return(input_list)
    
}
