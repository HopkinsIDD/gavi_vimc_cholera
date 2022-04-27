# Function updating proportion susceptible raster stack
# proportion susceptible is based on proportion vaccinated, considering both population migration (death) and vaccine waning effect (ve_direct)
update_sus_rasterStack <- function(datapath, 
                                   modelpath,
                                   country, 
                                   scenario,
                                   rawoutpath,
                                   pop,
                                   rc_list, 
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
