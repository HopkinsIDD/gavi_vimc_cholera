## Original incidence raster did not have a CRS so assigning one here and re-exporting it
## Jan 2021

## This was the original lambda raster stack that was used
# fn <- "../../svn/cholera-taxonomy/trunk/manuscripts/GAVI Impact Estimation/data/afro_2010-2016_lambdas_raster_stack_upd.grd"
# afr <- raster::stack(fn)
# raster::crs(afr) <- "+proj=longlat +datum=WGS84 +no_defs"

## exported 20x20 stack with correct CRS
# raster::writeRaster(afr, "input_data/incidence/afro_2010-2016_lambda.grd", overwrite=TRUE)

## exported mean 20x20 layer with correct CRS
# afr_mean <- raster::calc(afr, fun=mean, na.rm=TRUE)
# raster::writeRaster(afr_mean, "input_data/incidence/afro_2010-2016_lambda_mean.grd", overwrite=TRUE)

## finally decided we needed 1x1km raster stack to ensure best alignment
afr <- raster::stack("input_data/incidence/afro_2010-2016_lambda.grd")

afr5k <- raster::disaggregate(afr, 4, filename = "input_data/incidence/afro_2010-2016_lambda_5k.tif")

afr5k_mean <- raster::calc(afr5k, fun = mean, na.rm = TRUE, filename = "input_data/incidence/afro_2010-2016_lambda_5k_mean.tif")
