#' @name create_sus_modelInputs
#' @title create_sus_modelInputs
#' @description Create raster for proportion of population susceptible in each year, after accounting for vaccination
#' @param datapath path to input data 
#' @param modelpath path to montagu files
#' @param country country code
#' @param scenario Unique string that identifies the coverage scenario name
#' @param rawoutpath path to raw model output files
#' @param vacc_alloc object returned from [`allocate_vaccine()`]
#' @param ve_direct vaccine effect function (default: generate_pct_protect_function())
#' @param clean logical that indicates whether existing sus_files should be deleted
#' @return NULl
#' @export
create_sus_modelInputs <- function(
  datapath, 
  modelpath, 
  country, 
  scenario,
  rawoutpath,
  vacc_alloc,
  ve_direct,
  clean){
  
  ## create template and inputs
  years_ls <- get_model_years(modelpath, country, vacc_alloc)
  model_years <- years_ls[["model_years"]]
  output_years <- years_ls[["output_years"]]
  first_output_year <- output_years[1]
  real_output_years <- output_years[-1]

  startpop_raster <- create_model_pop_raster(datapath, modelpath, country, first_output_year)
  raster1_template <- raster::calc(startpop_raster, fun = function(x){ifelse(x>0, x/x, 0)})

  #### TRACK AND MODEL VACCINE EFFECTS AT GRID LEVEL ####
  ## create raster template for proportion susceptible 
  ## assume full susceptibility at the start of the model 
  dir.create(file.path(rawoutpath, scenario), showWarnings = FALSE)
  sus_out_fn <- paste0(rawoutpath, "/", scenario, "/", country, "_sus.tif")
  vacc_out_fn <- paste0(rawoutpath, "/", scenario, "/", country, "_vacc.tif")
  pop_out_fn <- paste0(rawoutpath, "/", scenario, "/", country, "_pop.tif")
  vacc_rasterStack <- raster::stack(vacc_out_fn)
  pop_rasterStack <- raster::stack(pop_out_fn)

  if (clean & file.exists(sus_out_fn)){
    message(paste("Clean existing", sus_out_fn))
    file.remove(sus_out_fn)
  }

  if (!file.exists(sus_out_fn)){

    ## create dummy to start stack (will be deleted later)
    sus_rasterStack <- raster::stack(raster1_template)
    
    ## add a new proportion suspectible layer to sus_rasterStack for each real model year
    for (j in 1:length(output_years)){ ## loop through model years

      ## Get life expectancy data from file
      mu <- 1/(import_country_lifeExpectancy_1yr(mpathname, country, output_years[j]))
      message(paste("Modeling susceptibility:", country, output_years[j], 1/mu, "life expectancy"))

      if (!(output_years[j] %in% model_years)){
        
        ## full susceptibility assumed before vaccination start and >8 years after vaccination ends (see get_model_years function)
        tmp <- raster1_template
        sus_rasterStack <- raster::addLayer(sus_rasterStack, tmp)

      } else{

        ## calculate remaining susceptible pop after vacc starts (loss due to migration and waning VE)
        tmp <- raster1_template
        popj <- raster::subset(pop_rasterStack, subset = j, drop = FALSE)

        for (k in 1:j){ ## loop through pre-j years
          popk <- raster::subset(pop_rasterStack, subset = k, drop = FALSE)
          vacck <- raster::subset(vacc_rasterStack, subset = k, drop = FALSE)

          pkj <- raster::overlay(popk, popj,
            fun = function(x, y){x*(1-((j-k)*mu))/y}
            ) ## population retention (measures turnover rate due to death) from years k into year j
          ve_j_k <- as.numeric(ve_direct(j-k+1)) ## couldn't put the ve_direct function in raster::overlay even though it just returns a scalar
          prob_still_protected <- raster::overlay(vacck, pkj, 
            fun = function(x, y){return(x*y*ve_j_k)} 
            ) 
          tmp <- raster::overlay(tmp, prob_still_protected,
            fun = function(x, y){x*(1-y)}
            )

          ## track proportion of full pop immune for cases averted calculation
          prop_protected <- raster::cellStats(
            raster::overlay(prob_still_protected, popj,
              fun = function(x, y){x*y}
              ), "sum")/raster::cellStats(popj, "sum") 

          track_prop_immune <- c(country = country,
                                scenario = scenario,
                                year = real_output_years[j],
                                campaignYear = real_output_years[k],
                                prop_immune_in_year = prop_protected)
          # print(dplyr::glimpse(track_prop_immune))

          rm(popk, vacck, pkj, prob_still_protected)
          gc()
        } #endkfor

        sus_rasterStack <- raster::addLayer(sus_rasterStack, tmp)
        rm(tmp, popj)
        gc()

      }
   
    } #endjfor

    ## drop initial sus rast layer, which was a dummy
    sus_rasterStack <- raster::dropLayer(sus_rasterStack, 1)
    if (dim(sus_rasterStack)[3] != length(output_years)){
      stop(paste("susceptibility rasterStack has", dim(sus_rasterStack)[3], "layers but there are", length(output_years), "model years."))
    }

    message(paste("Write", sus_out_fn))
    sus_rs <- align_rasters(datapath, country, sus_rasterStack)
    raster::writeRaster(sus_rs, filename = sus_out_fn)

    rm(startpop_raster, raster1_template, sus_rasterStack, sus_rs)
    gc()

  } else {
    message(paste("Skip creation", sus_out_fn))
  }

  return(NULL)
}
  