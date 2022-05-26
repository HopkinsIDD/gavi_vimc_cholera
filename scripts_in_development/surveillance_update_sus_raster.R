# Function updating proportion susceptible raster stack
# proportion susceptible is based on proportion vaccinated, considering both population migration (death) and vaccine waning effect (ve_direct)
update_sus_rasterStack_optimized <- function( datapath, 
                                              modelpath,
                                              country, 
                                              scenario,
                                              rawoutpath,
                                              pop, 
                                              ve_direct,
                                              baseline_year,
                                              model_year,
                                              last_vac_ef_year = 5   # empty list
                                              ){
  message(paste0("Stacking proportion susceptible raster layer of ", model_year, "."))
  
  ### First make sure that the last vaccine effective year is the 5th year 
  if(! (ve_direct(last_vac_ef_year)>0 & ve_direct(last_vac_ef_year+0.001)[1]==0) ){
    stop("The last vaccine effective year is not the 5th year, please check the update_sus_rasterStack_optimized function and change the preset value. ")
  }
  

  # ### Create the file names 
  # vac_admin1_fn <- paste0("intermediate_raster/", country, "_vac_admin1_", model_year, ".tif")
  # vac_admin2_fn <- paste0("intermediate_raster/", country, "_vac_admin2_", model_year, ".tif")
  # vac_pop_fn <- paste0("intermediate_raster/", country, "_vac_pop_", model_year, ".tif")

  
  # vac_pop_fn <- paste0("intermediate_raster/", country, "_vac_pop_", model_year, ".tif")
  # pop_rasterStack <- raster::stack(vac_pop_fn)
  raster1_template <- raster::calc(pop, fun = function(x){ifelse(!is.na(x), 1, NA)})
  j <- model_year - baseline_year + 1
  vac_pop_fn <- paste0("intermediate_raster/", country, "_vac_pop_", model_year, ".tif")
  popj <- raster::stack(vac_pop_fn) 
  tmp1 <- raster1_template
  tmp2 <- raster1_template

  # get life expectancy data from file       
  mu <- 1/(import_country_lifeExpectancy_1yr(modelpath, country, model_year))
  message(paste("Modeling susceptibility:", country, model_year, 1/mu, "life expectancy"))

  for (year in max(baseline_year, model_year-last_vac_ef_year+1):model_year){
    
    k <- year - baseline_year + 1 
    vac_admin1_fn <- paste0("intermediate_raster/", country, "_vac_admin1_", year, ".tif")
    vac_admin2_fn <- paste0("intermediate_raster/", country, "_vac_admin2_", year, ".tif")
    vac_pop_fn <- paste0("intermediate_raster/", country, "_vac_pop_", year, ".tif")
    popk <- raster::stack(vac_pop_fn)
    vacck_admin1 <- raster::stack(vac_admin1_fn)
    vacck_admin2 <- raster::stack(vac_admin2_fn)
      
    # population retention (measures turnover rate due to death) from years k into year j -- this assumes that the pop during one same year is constant**********
    # a fix for the 0's on both population rasters (0 -> 1), good thing is that the changed population rasters will not be used elsewhere 
    popk[popk == 0] <- 1
    popj[popj == 0] <- 1
    pkj <- raster::overlay(popk, popj, fun = function(x, y){x*(1-((j-k)*mu))/y}) 
    ve_j_k <- as.numeric(ve_direct(j-k+1))
      
    prob_still_protected_admin1 <- raster::overlay(vacck_admin1, pkj, fun = function(x, y){return(x*y*ve_j_k)}) 
    prob_still_protected_admin2 <- raster::overlay(vacck_admin2, pkj, fun = function(x, y){return(x*y*ve_j_k)}) 
      
    # get the new sus raster layer -- from protected to still susceptible 
    tmp1 <- raster::stack(raster::overlay(tmp1, prob_still_protected_admin1, fun = function(x, y){x*(1-y)}))
    tmp2 <- raster::stack(raster::overlay(tmp2, prob_still_protected_admin2, fun = function(x, y){x*(1-y)}))
      
    rm(popk, popj, vacck_admin1, vacck_admin2, pkj, prob_still_protected_admin1, prob_still_protected_admin2)
    gc()
      
  }
    
  sus_list <- list( "sus_rasterStack_admin1" = tmp1,
                    "sus_rasterStack_admin2" = tmp2)
    
  return(sus_list)

}



save_sus_raster <- function(datapath, modelpath, country, nsamples, model_year, sus_list, 
                            rawoutpath = NULL, 
                            clean = NULL, 
                            scenario = NULL, 
                            incidence_rate_trend = NULL, 
                            outbreak_multiplier = NULL, 
                            vac_incid_threshold = NULL, 
                            surveillance_scenario = NULL){

  ### Create the file names 
  sus_admin1_fn <- paste0("intermediate_raster/", country, "_sus_admin1_", model_year, ".tif")
  sus_admin2_fn <- paste0("intermediate_raster/", country, "_sus_admin2_", model_year, ".tif")
  
  sus1_out_fn <- paste0(rawoutpath, "/", scenario, "/", paste("incid", incidence_rate_trend, "outbk", outbreak_multiplier, 
                        vac_incid_threshold, surveillance_scenario, country, sep = "_"), "_sus_admin1_", model_year, ".tif")
  sus2_out_fn <- paste0(rawoutpath, "/", scenario, "/", paste("incid", incidence_rate_trend, "outbk", outbreak_multiplier, 
                        vac_incid_threshold, surveillance_scenario, country, sep = "_"), "_sus_admin2_", model_year, ".tif")
  

  ### Sometimes convert the structure 
  if(!is.null(sus_list) & is.null(names(sus_list)) & !is.null(rawoutpath)){
    sus_admin1 <- sus_list[[1]]$sus_rasterStack_admin1
    sus_admin2 <- sus_list[[1]]$sus_rasterStack_admin2
    
    for(layer_idx in 2:nsamples){
      sus_admin1 <- raster::stack(sus_admin1, sus_list[[layer_idx]]$sus_rasterStack_admin1)
      sus_admin2 <- raster::stack(sus_admin2, sus_list[[layer_idx]]$sus_rasterStack_admin2)
    }

    message(paste("Writing proportion susceptible rasterStack for", country))
    dir.create(paste0(rawoutpath, "/", scenario, "/"), showWarnings = FALSE)
    if( !file.exists(sus1_out_fn) | (file.exists(sus1_out_fn)&clean) ){
      raster::writeRaster(sus_admin1, filename = sus1_out_fn, overwrite = TRUE)}
    if( !file.exists(sus2_out_fn) | (file.exists(sus2_out_fn)&clean) ){
      raster::writeRaster(sus_admin2, filename = sus2_out_fn, overwrite = TRUE)}

    return(NULL)
  }


  ### Save 
  if(is.null(rawoutpath)){
    dir.create(paste0("intermediate_raster/"), showWarnings = FALSE)
    raster::writeRaster(sus_list$sus_rasterStack_admin1, filename = sus_admin1_fn, overwrite = TRUE)
    raster::writeRaster(sus_list$sus_rasterStack_admin2, filename = sus_admin2_fn, overwrite = TRUE)
  }else{
    message(paste("Writing proportion susceptible rasterStack for", country))
    dir.create(paste0(rawoutpath, "/", scenario, "/"), showWarnings = FALSE)
    if( !file.exists(sus1_out_fn) | (file.exists(sus1_out_fn)&clean) ){
      file.rename(from=sus_admin1_fn, to=sus1_out_fn)}
    if( !file.exists(sus2_out_fn) | (file.exists(sus2_out_fn)&clean) ){
      file.rename(from=sus_admin2_fn, to=sus2_out_fn)}
  }

}



update_sus_rasterStack <- function(datapath, 
                                   modelpath,
                                   country, 
                                   scenario,
                                   rawoutpath,
                                   pop, 
                                   ve_direct,
                                   baseline_year,
                                   model_year, 
                                   input_list, # generated from update_input_rasterStack() function
                                   sus_list   # rasterStack of proportion susceptible generated in last year
                                  ){
  
  ### Get the rasters ready 
  pop_rasterStack <- input_list[["pop_rasterStack"]]
  vacc_rasterStack_admin1 <- input_list[["vacc_rasterStack_admin1"]]
  vacc_rasterStack_admin2 <- input_list[["vacc_rasterStack_admin2"]]

  raster1_template <- raster::calc(pop, fun = function(x){ifelse(!is.na(x), 1, NA)})
    
  message(paste0("Stacking proportion susceptible raster layer of ", model_year, "."))
  
  if(!is.null(sus_list)){
    sus_rasterStack_admin1 <- sus_list[["sus_rasterStack_admin1"]]
    sus_rasterStack_admin2 <- sus_list[["sus_rasterStack_admin2"]]
  }
  
    
  j = model_year - baseline_year + 1 #the j'th year into simulation 
  popj <- raster::subset(pop_rasterStack, subset = j, drop = FALSE)
  tmp1 <- raster1_template
  tmp2 <- raster1_template
  
  # get life expectancy data from file       
  mu <- 1/(import_country_lifeExpectancy_1yr(modelpath, country, model_year))
  message(paste("Modeling susceptibility:", country, model_year, 1/mu, "life expectancy"))

  ### Loop through the years before the current running year (including the current model year)
  for(k in 1:j){ #may not need that many years given the protection of vaccines is not to last for more than 5 years 
      
    popk <- raster::subset(pop_rasterStack, subset = k, drop = FALSE)
    vacck_admin1 <- raster::subset(vacc_rasterStack_admin1, subset = k, drop = FALSE)
    vacck_admin2 <- raster::subset(vacc_rasterStack_admin2, subset = k, drop = FALSE)
      
    # population retention (measures turnover rate due to death) from years k into year j -- this assumes that the pop during one same year is constant**********
    # a fix for the 0's on both population rasters (0 -> 1), good thing is that the changed population rasters will not be used elsewhere 
    popk[popk == 0] <- 1
    popj[popj == 0] <- 1
    pkj <- raster::overlay(popk, popj, fun = function(x, y){x*(1-((j-k)*mu))/y}) 
    # pkj <- popk*(1-((j-k)*mu))/popj #new
    ve_j_k <- as.numeric(ve_direct(j-k+1))
      
    prob_still_protected_admin1 <- raster::overlay(vacck_admin1, pkj, fun = function(x, y){return(x*y*ve_j_k)}) 
    prob_still_protected_admin2 <- raster::overlay(vacck_admin2, pkj, fun = function(x, y){return(x*y*ve_j_k)}) 
      
    # get the new sus raster layer -- from protected to still susceptible 
    tmp1 <- raster::overlay(tmp1, prob_still_protected_admin1, fun = function(x, y){x*(1-y)})
    tmp2 <- raster::overlay(tmp2, prob_still_protected_admin2, fun = function(x, y){x*(1-y)})
      
    rm(popk, vacck_admin1, vacck_admin2, pkj, prob_still_protected_admin1, prob_still_protected_admin2)
    gc()
      
  } #endkfor
    
  ### Append and store and return the rasters 
  # append new sus raster layers
  if(is.null(sus_list)){
    sus_rasterStack_admin1 <- tmp1
    sus_rasterStack_admin2 <- tmp2
  }else{
    sus_rasterStack_admin1 <- raster::addLayer(sus_rasterStack_admin1, tmp1)
    sus_rasterStack_admin2 <- raster::addLayer(sus_rasterStack_admin2, tmp2)
  }
  rm(tmp1, tmp2, popj)
  gc()
  
  # align raster to the extent and resolution of the worldpop rasters and GADM shapefiles -- necessary? 
  sus_rs_admin1 <- sus_rasterStack_admin1
  sus_rs_admin2 <- sus_rasterStack_admin2
  # sus_rs_admin1 <- align_rasters(datapath, country, sus_rasterStack_admin1)
  # sus_rs_admin2 <- align_rasters(datapath, country, sus_rasterStack_admin2)
    
  rm(raster1_template, pop_rasterStack, vacc_rasterStack_admin1, vacc_rasterStack_admin2,
      sus_rasterStack_admin1, sus_rasterStack_admin2)
  gc()
    
  sus_list <- list( "sus_rasterStack_admin1" = sus_rs_admin1,
                    "sus_rasterStack_admin2" = sus_rs_admin2)
    
  return(sus_list)

}



### The following function will be called if the assumption that the vacc begins and finishes at the start of the year and affects the whole year's cases 
# function to load the susceptible rasterStack list for the first year
create_first_year_sus <- function(datapath, modelpath, country, scenario, rawoutpath,
                                  rc_list, ve_direct,
                                  baseline_year, 
                                  pop,
                                  input_list, # generated from update_input_rasterStack() function
                                  clean){
    
    raster1_template <- raster::calc(pop, fun = function(x){ifelse(x>0, x/x, 0)})
    message("Loading the first layer of proportion susceptible rasterStack.")

    # assume full susceptibility in the first year
    sus_rasterStack <- raster::stack(raster1_template)
    
    sus_list <- list("sus_rasterStack_admin1" = sus_rasterStack,
                     "sus_rasterStack_admin2" = sus_rasterStack)
 
    return(sus_list)
    rm(raster1_template, sus_rasterStack)
}
