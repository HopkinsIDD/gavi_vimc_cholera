require(devtools)
install_version("raster", version = "3.4.13", repos = "http://cran.us.r-project.org")
install_version("exactextractr", version = "0.6.1", repos = "http://cran.us.r-project.org")

datapath <- "/home/kaiyuezou/VIMC_Model/surveillance_project/gavi_vimc_cholera/input_data"
country <- "ETH"
afr <- raster::raster(paste0(datapath, "/incidence/afro_2010-2016_lambda_5k_mean.tif"))
shp1 <- GADMTools::gadm_sf_loadCountries(c(country), level = 1, basefile = file.path(datapath, "shapefiles/"))$sf

incid1 <- exactextractr::exact_extract(afr, shp1, 'mean')
