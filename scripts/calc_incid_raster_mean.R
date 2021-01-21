## Original incidence raster did not have a CRS so assigning one here and re-exporting it
## Jan 2021
fn <- "C:/Users/eclee/Dropbox (Bansal Lab)/IDDynamicsJHU/svn/cholera-taxonomy/trunk/manuscripts/GAVI Impact Estimation/data/afro_2010-2016_lambdas_raster_stack_upd.grd"
afr <- raster::stack(fn)
raster::crs(afr) <- "+proj=longlat +datum=WGS84 +no_defs"

raster::writeRaster(afr, "C:/Users/eclee/Dropbox (Bansal Lab)/IDDynamicsJHU/gh/gavi_vimc_cholera/input_data/incidence/afro_2010-2016_lambda.grd", overwrite=TRUE)
  
afr_mean <- raster::calc(afr, fun=mean, na.rm=TRUE)
# raster::writeRaster(afr_mean, "C:/Users/eclee/Dropbox (Bansal Lab)/IDDynamicsJHU/gh/gavi_vimc_cholera/input_data/incidence/afro_2010-2016_lambda_mean.grd", overwrite=TRUE)
