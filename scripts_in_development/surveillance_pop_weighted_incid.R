pop_weighted_admin_mean_incid <- function(datapath, modelpath, incidence_rate_raster, pop_raster, country, admin_shp){

  # afr <- raster::raster(paste0(datapath, "/incidence/afro_2010-2016_lambda_5k_mean.tif"))
  shp0 <- load_shp0_by_country(datapath, country)
  country_incid_raster <- raster::mask(raster::crop(raster::stack(incidence_rate_raster), shp0, snap = "out"), shp0, updatevalue = NA)
  # pop_certain_year <- ocvImpact::create_model_pop_raster(datapath, modelpath, country, year)

  pop_for_weights <- raster::resample(pop_raster, country_incid_raster)
  # pop_for_weights <- pop_certain_year
  incid <- exactextractr::exact_extract(incidence_rate_raster, admin_shp, 'weighted_mean', weights = area(pop_for_weights))
  # incid <- exactextractr::exact_extract(country_incid_raster, admin_shp, "mean")
  # incid <- exactextractr::exact_extract(country_incid_raster, admin_shp, 'weighted_mean', weights = pop_for_weights)

  return(incid)
}
