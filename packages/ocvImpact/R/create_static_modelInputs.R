#' @name create_static_modelInputs
#' @title create_static_modelInputs
#' @description Create proportion of population vaccinated and total population rasterStacks. Write vaccination raster to file and export population raster
#' @param datapath path to input data
#' @param modelpath path to montagu files
#' @param country country code
#' @param scenario Unique string that identifies the coverage scenario name
#' @param rawoutpath path to raw model output files
#' @param vacc_alloc object returned from [`allocate_vaccine()`]
#' @param cache montagu cache
#' @param clean logical that indicates whether existing model output files (vacc & pop) should be deleted
#' @return rasterStack of total population for model years
#' @export
#' @include utils.R
create_static_modelInputs <- function(
  datapath,
  modelpath,
  country,
  scenario,
  rawoutpath,
  vacc_alloc,
  cache,
  clean){

  ## create template and inputs
  years_ls <- get_model_years(modelpath, country, vacc_alloc, cache)
  output_years <- years_ls[["output_years"]]
  first_output_year <- output_years[1]

  startpop_raster <- create_model_pop_raster(datapath, modelpath, country, first_output_year, cache)
  raster0_template <- raster::calc(startpop_raster, fun = function(x){x*0})

  #### CREATE proportion vaccinated AND total population RASTERS FOR MODELING ####

  ## create rasterStack for
  ## 1. proportion of population vaccinated in each grid cell
  ## 2. population
  vacc_rasterStack <- raster::stack(raster0_template)
  pop_rasterStack <- raster::stack(startpop_raster)
  dir.create(file.path(rawoutpath, scenario), showWarnings = FALSE)
  vacc_out_fn <- paste0(rawoutpath, "/", scenario, "/", country, "_vacc.tif")
  pop_out_fn <- paste0(rawoutpath, "/", scenario, "/", country, "_pop.tif")

  if(clean & file.exists(vacc_out_fn)){
    message(paste("Clean existing", vacc_out_fn))
    file.remove(vacc_out_fn)
  }
  if(clean & file.exists(pop_out_fn)){
    message(paste("Clean existing", pop_out_fn))
    file.remove(pop_out_fn)
  }

  if (file.exists(vacc_out_fn) & file.exists(pop_out_fn)){ ## if vacc raster file already exists for this scenario
    message(paste("Skip creation", vacc_out_fn, pop_out_fn))

  } else { ## write vacc & pop raster to file for this scenario

    if(file.exists(vacc_out_fn)){
      message(paste("Clean existing", vacc_out_fn))
      file.remove(vacc_out_fn)
    }
    if(file.exists(pop_out_fn)){
      message(paste("Clean existing", pop_out_fn))
      file.remove(pop_out_fn)
    }

    for (year in output_years){
      message(paste("In", year, "of static_modelInputs"))

      if (year %in% unique(vacc_alloc$vacc_year)){
        new_layer <- fasterize::fasterize(
          dplyr::filter(vacc_alloc, vacc_year == year),
          raster0_template,
          field = "actual_prop_atleast_1dose_vaccinated",
          fun = "last",
          background = 0
          )
        new_vacc_layer <- raster::mask(new_layer, raster0_template, updatevalue = NA)
      } else{
        ## add 0 layer if there was no vaccination that year
        new_vacc_layer <- raster0_template
      } #endifelse

      vacc_rasterStack <- raster::addLayer(vacc_rasterStack, new_vacc_layer)

      new_pop_layer <- create_model_pop_raster(datapath, modelpath, country, year, cache)
      pop_rasterStack <- raster::addLayer(pop_rasterStack, new_pop_layer)

      rm(new_vacc_layer, new_pop_layer)
      gc()

    } #endfor model_years

    ## drop initial vacc & pop rast layers, which were dummies
    vacc_rs <- raster::dropLayer(vacc_rasterStack, 1)
    pop_rs <- raster::dropLayer(pop_rasterStack, 1)

    if ((dim(vacc_rs)[3] != length(output_years)) &
        (dim(pop_rs)[3] != length(output_years))){

      stop(paste(
        "The vaccine intervention is implemented in",
        dim(vacc_rs)[3],
        "years, the population data is stacked for",
        dim(pop_rs)[3],
        "years, but the model will run for",
        length(output_years),
        "years."))
    }

    # vacc_rs <- align_rasters(datapath, country, vacc_rasterStack) ## may be obsolete 1/21
    message(paste("Write", vacc_out_fn))
    raster::writeRaster(vacc_rs, filename = vacc_out_fn)

    # pop_rs <- align_rasters(datapath, country, pop_rasterStack) ## may be obsolete 1/21
    message(paste("Write", pop_out_fn))
    raster::writeRaster(pop_rs, filename = pop_out_fn)
  }

  rm(vacc_rasterStack, startpop_raster, raster0_template, pop_rasterStack)
  gc()

  return(NULL)

}
