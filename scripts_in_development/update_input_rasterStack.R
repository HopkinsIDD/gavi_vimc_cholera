# Function: update_input_rasterStack() to get the input rasterstack
update_input_rasterStack <- function(datapath,
                               modelpath, 
                               country, 
                               scenario,
                               rawoutpath,
                               rc_list, 
                               baseline_year,
                               this_year,
                               pop, # pop raster for the latest year
                               input_list # input an empty list for the first year, for the following year, input is that list from last year
                               ){

    
  message(paste("Stacking population and vaccination proportion raster layer of", this_year, "."))
    
  # read in stacks calculated from last campaign year
  vacc_rasterStack_admin1 <- input_list[["vacc_rasterStack_admin1"]]
  vacc_rasterStack_admin2 <- input_list[["vacc_rasterStack_admin2"]]
  pop_rasterStack <- input_list[["pop_rasterStack"]]
  raster0_template <- raster::calc(pop, fun = function(x){x*0})
  
    
  # append new layer
  pop_rasterStack <- raster::addLayer(pop_rasterStack, pop)

  shp1_targeted <- rc_list[[1]] %>% filter(target_year == this_year & is_target == 1)
  shp2_targeted <- rc_list[[2]] %>% filter(target_year == this_year & is_target == 1)
  
  # check point: if there is not a single place targeted this year.
  if(nrow(shp1_targeted) == 0){
    
    message(paste("There is no admin 1 level areas vaccinated in", this_year, "in", country, "."))
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
    
    message(paste("There is no admin 2 level areas vaccinated in", this_year, "in", country, "."))
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
  
  vacc_rasterStack_admin1 <- raster::addLayer(vacc_rasterStack_admin1, new_vacc_layer_admin1) 
  vacc_rasterStack_admin2 <- raster::addLayer(vacc_rasterStack_admin2, new_vacc_layer_admin2)

  input_list = list("vacc_rasterStack_admin1" = vacc_rasterStack_admin1,
                    "vacc_rasterStack_admin2" = vacc_rasterStack_admin2,
                    "pop_rasterStack" = pop_rasterStack)
  rm(new_vacc_layer_admin1, new_vacc_layer_admin2)

  return(input_list)
  
}



create_first_year_input <- function(datapath, modelpath, country,
                                  baseline_year, rc_list,
                                  pop){

    message("Loading first layer of population and vaccination proportion rasterStack.")
    
    raster0_template <- raster::calc(pop, fun = function(x){x*0})
    shp1_targeted <- rc_list[[1]] %>% filter(target_year == baseline_year & is_target == 1)
    shp2_targeted <- rc_list[[2]] %>% filter(target_year == baseline_year & is_target == 1)

    # check point: if there is not a single place targeted this year.
    if(nrow(shp1_targeted) == 0){
      
      message(paste("There is no admin 1 level areas vaccinated in", this_year, "in", country, "."))
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
      
      message(paste("There is no admin 2 level areas vaccinated in", this_year, "in", country, "."))
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