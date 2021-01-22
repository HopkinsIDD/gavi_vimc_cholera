#' @name create_expectedCases
#'@title create_expectedCases
#' @description Create proportion of population vaccinated and total population rasterStacks. Write vaccination raster to file and export population raster
#' @param datapath path to input data 
#' @param modelpath path to montagu files
#' @param country country code
#' @param scenario Unique string that identifies the coverage scenario name
#' @param rawoutpath path to raw model output files
#' @param vacc_alloc object returned from [`allocate_vaccine()`]
#' @param indirect_mult function for calculating indirect effects
#' @param secular_trend_mult function for calculating secular trends in mean annual incidence
#' @param nsamples number of stochastic samples to use
#' @param is_cf logical indicating whether to run counterfactual model
#' @param clean logical indicating whether to delete existing files; pass to [`create_incid_raster()`]
#' @return dataframe with 
#' @export
create_expectedCases <- function(
  datapath,
  modelpath,
  country,
  scenario,
  rawoutpath,
  vacc_alloc,
  indirect_mult,
  secular_trend_mult,
  nsamples,
  is_cf,
  clean
  ){

  ## create template and inputs
  years_ls <- get_model_years(modelpath, country, vacc_alloc)
  model_years <- years_ls[["model_years"]]
  if (is_cf){
    model_years <- NULL
  }
  tmp <- import_centralburden_template(mpathname, country)
  output_years <- sort(unique(tmp$year))

  sus_out_fn <- paste0(rawoutpath, "/", scenario, "/", country, "_sus.tif")
  vacc_out_fn <- paste0(rawoutpath, "/", scenario, "/", country, "_vacc.tif")
  pop_out_fn <- paste0(rawoutpath, "/", scenario, "/", country, "_pop.tif")

  ## write to file and import cholera incidence estimates
  lambda <- create_incid_raster(datapath, country, rawoutpath, nsamples, clean)
  sus_rasterStack <- raster::brick(sus_out_fn)
  vacc_rasterStack <- raster::brick(vacc_out_fn)
  pop_rasterStack <- raster::brick(pop_out_fn)
  shp0 <- load_shapefile_by_country(datapath, country, simple = TRUE)

  ec_ls <- lapply(output_years, function(oy){

    yr_index <- which(oy == output_years)
    pop_rasterLayer <- raster::subset(pop_rasterStack, yr_index, drop = FALSE)
    vacc_rasterLayer <- raster::subset(vacc_rasterStack, yr_index, drop = FALSE)
    sus_rasterLayer <- raster::subset(sus_rasterStack, yr_index, drop = FALSE)

    if (oy %in% model_years){ ## consecutive years where vaccine dynamics are in play

      ## make new indirect effects template
      indirect_rasterLayer <- pop_rasterLayer
      raster::values(indirect_rasterLayer) <- indirect_mult(1-as.numeric(raster::values(sus_rasterLayer)))

      ec_rasterStack <- raster::overlay(
        sus_rasterLayer,
        pop_rasterLayer,
        lambda,
        indirect_rasterLayer,
        fun = function(x, y, z, a){
          x*y*z*a*secular_trend_mult(year = oy)
        },
        recycle = TRUE)

    } else{
      
      ec_rasterStack <- raster::overlay(
        pop_rasterLayer,
        lambda,
        fun = function(x, y){
          return(x*y*secular_trend_mult(year = oy))
        },
        recycle = TRUE, unstack = TRUE)

    }

    ec_yr <- exactextractr::exact_extract(ec_rasterStack, shp0, fun = "sum", stack_apply = TRUE)
    mean_incid <- exactextractr::exact_extract(lambda, shp0, fun = "weighted_mean", weights = pop_rasterLayer, stack_apply = TRUE)

    ec_vec <- as.numeric(ec_yr)
    mean_incid_vec <- as.numeric(mean_incid)
    rc <- tibble::tibble(country = country, year = oy, run_id = seq_along(ec_vec), ec = ec_vec, incid_rate = mean_incid_vec)

    return(rc)

  })

  ec_final <- data.table::rbindlist(ec_ls)
  
  rm(lambda, sus_rasterStack, vacc_rasterStack, pop_rasterStack, shp, ec_ls)
  gc()

  return(ec_final)
}