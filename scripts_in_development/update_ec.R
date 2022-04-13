# Function to update the expected cases rasterStack
update_ec <- function(datapath,
                       modelpath,
                       country,
                       scenario,
                       rawoutpath,
                       rc_list, 
                       baseline_year,
                       sus_list, # susceptible rasterStack generated in the last function
                       input_list,      # population and proportion vaccinated rasterStack generated in this loop
                       indirect_mult,
                       secular_trend_mult,
                       lambda, # loaded baseline incidence of this country, using: lambda <- create_incid_raster(modelpath, datapath, country, nsamples, redraw)
                       ec_list # list of expected cases raster from last loop
                       # nsamples, 
                       # is_cf,
                       # redraw
                      ){
  
  this_year <- max(rc_list[[1]]$target_year)
  # tmp <- import_centralburden_template(modelpath, country, redownload = FALSE)
  
  message(paste0("Stacking expected cases raster layer of ", this_year, "."))

  #pop_rasterStack <- raster::brick(input_list[["pop_rasterStack"]])
  #vacc_rasterStack_admin1 <- raster::brick(input_list[["vacc_rasterStack_admin1"]])
  #vacc_rasterStack_admin2 <- raster::brick(input_list[["vacc_rasterStack_admin2"]])
  #sus_rasterStack_admin1 <- raster::brick(sus_list[["sus_rasterStack_admin1"]])
  #sus_rasterStack_admin2 <- raster::brick(sus_list[["sus_rasterStack_admin2"]])

  # subset raster layers of this_year
  yr_index <- which(this_year == unique(rc_list[[1]]$target_year))
  pop_rasterLayer <- raster::subset(input_list[["pop_rasterStack"]], yr_index, drop = FALSE)
  vacc_rasterLayer_admin1 <- raster::subset(input_list[["vacc_rasterStack_admin1"]], yr_index, drop = FALSE)
  vacc_rasterLayer_admin2 <- raster::subset(input_list[["vacc_rasterStack_admin2"]], yr_index, drop = FALSE)
  sus_rasterLayer_admin1 <- raster::subset(sus_list[["sus_rasterStack_admin1"]], yr_index, drop = FALSE)
  sus_rasterLayer_admin2 <- raster::subset(sus_list[["sus_rasterStack_admin2"]], yr_index, drop = FALSE)
  
  ## make new indirect effects template
  indirect_rasterLayer_admin1 <- pop_rasterLayer
  indirect_rasterLayer_admin2 <- pop_rasterLayer
  raster::values(indirect_rasterLayer_admin1) <- indirect_mult(1-as.numeric(raster::values(sus_rasterLayer_admin1)))
  raster::values(indirect_rasterLayer_admin2) <- indirect_mult(1-as.numeric(raster::values(sus_rasterLayer_admin2)))
  incid_trend <- secular_trend_mult(year = this_year)
    
  # calculate ec at admin1 (the output ec_admin1 has 30 layers, because nsample = 30)
  ec_admin1 <- 
      raster::overlay(
        sus_rasterLayer_admin1,
        pop_rasterLayer,
        lambda,
        indirect_rasterLayer_admin1,
        fun = function(x, y, z, a){
          x*y*z*a*incid_trend
        })

  ec_admin2 <- 
    raster::overlay(
        sus_rasterLayer_admin2,
        pop_rasterLayer,
        lambda,
        indirect_rasterLayer_admin2,
        fun = function(x, y, z, a){
          x*y*z*a*incid_trend
        })
    
  # calculate mean ec across the 30 layers or ec (other ways to do this?)
  ec_rasterLayer_admin1 <- raster::calc(ec_admin1, fun = mean, na.rm = T)
  ec_rasterLayer_admin2 <- raster::calc(ec_admin2, fun = mean, na.rm = T)
      
  # append new ec raster layer
  ec_list[["ec_rasterStack_admin1"]] <- raster::addLayer(ec_list[["ec_rasterStack_admin1"]], ec_rasterLayer_admin1)
  ec_list[["ec_rasterStack_admin2"]] <- raster::addLayer(ec_list[["ec_rasterStack_admin2"]], ec_rasterLayer_admin2)
  
  rm(ec_rasterLayer_admin1, ec_rasterLayer_admin2, 
     pop_rasterLayer, vacc_rasterLayer_admin1, vacc_rasterLayer_admin2,
     sus_rasterLayer_admin1, sus_rasterLayer_admin2)
  gc()
  
  return(ec_list)
  
}


create_first_year_ec <- function(datapath, modelpath, country, scenario, rawoutpath,
                                 rc_list, # incidence sf
                                 baseline_year,
                                 sus_list, # susceptible rasterStack generated in the last function
                                 input_list,      # population and proportion vaccinated rasterStack generated in this loop
                                 indirect_mult,
                                 secular_trend_mult,
                                 lambda # loaded baseline incidence of this country, using: lambda <- create_incid_raster(modelpath, datapath, country, nsamples, redraw)
                                 # nsamples, 
                                 # is_cf,
                                 # redraw
                                 ){
  
  message("Loading the first layer of expected cases rasterStack.")

  #read in input rasterStacks
  #pop_rasterStack <- raster::brick(input_list[["pop_rasterStack"]])
  #vacc_rasterStack_admin1 <- raster::brick(input_list[["vacc_rasterStack_admin1"]])
  #vacc_rasterStack_admin2 <- raster::brick(input_list[["vacc_rasterStack_admin2"]])
  #sus_rasterStack_admin1 <- raster::brick(sus_list[["sus_rasterStack_admin1"]])
  #sus_rasterStack_admin2 <- raster::brick(sus_list[["sus_rasterStack_admin2"]])

  # subset raster layers of this_year
  pop_rasterLayer <- raster::subset(input_list[["pop_rasterStack"]], 1, drop = FALSE)
  vacc_rasterLayer_admin1 <- raster::subset(input_list[["vacc_rasterStack_admin1"]], 1, drop = FALSE)
  vacc_rasterLayer_admin2 <- raster::subset(input_list[["vacc_rasterStack_admin2"]], 1, drop = FALSE)
  sus_rasterLayer_admin1 <- raster::subset(sus_list[["sus_rasterStack_admin1"]], 1, drop = FALSE)
  sus_rasterLayer_admin2 <- raster::subset(sus_list[["sus_rasterStack_admin2"]], 1, drop = FALSE)
    
  ## make new indirect effects template
  indirect_rasterLayer_admin1 <- pop_rasterLayer
  indirect_rasterLayer_admin2 <- pop_rasterLayer
  raster::values(indirect_rasterLayer_admin1) <- indirect_mult(1-as.numeric(raster::values(sus_rasterLayer_admin1)))
  raster::values(indirect_rasterLayer_admin2) <- indirect_mult(1-as.numeric(raster::values(sus_rasterLayer_admin2)))
  incid_trend <- secular_trend_mult(year = baseline_year)
    
  # calculate ec at admin1 (the output ec_admin1 has 30 layers, because nsample = 30)
  ec_admin1 <- 
      raster::overlay(
        sus_rasterLayer_admin1,
        pop_rasterLayer,
        lambda,
        indirect_rasterLayer_admin1,
        fun = function(x, y, z, a){
          x*y*z*a*incid_trend
        })

  ec_admin2 <- 
      raster::overlay(
        sus_rasterLayer_admin2,
        pop_rasterLayer,
        lambda,
        indirect_rasterLayer_admin2,
        fun = function(x, y, z, a){
          x*y*z*a*incid_trend
        })
  
  # calculate mean ec across the 30 layers or ec
  ec_rasterLayer_admin1 <- raster::calc(ec_admin1, fun = mean, na.rm = T)
  ec_rasterLayer_admin2 <- raster::calc(ec_admin2, fun = mean, na.rm = T)
    
  # stack the first layer
  ec_rasterStack_admin1 <- raster::stack(ec_rasterLayer_admin1)
  ec_rasterStack_admin2 <- raster::stack(ec_rasterLayer_admin2)

  ec_list <- list("ec_rasterStack_admin1" = ec_rasterStack_admin1,
                  "ec_rasterStack_admin2" = ec_rasterStack_admin2)
  
  rm(ec_rasterLayer_admin1, ec_rasterLayer_admin2, 
     pop_rasterLayer, vacc_rasterLayer_admin1, vacc_rasterLayer_admin2,
     sus_rasterLayer_admin1, sus_rasterLayer_admin2)
  gc()
  
  return(ec_list)

}