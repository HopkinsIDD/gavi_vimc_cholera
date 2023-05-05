test_that("create_incid_raster works", {
  
  # Try getting the root directory for the convenience of loading files 
  home_dir <- gsub("packages/ocvImpact/tests/testthat", "", getwd())

  # Regular African countries with raster data work 
  library(dplyr)
  tmpfile <- tempfile(fileext = ".yml")
  yaml::write_yaml(list(setting = list(random_seed = 1)), tmpfile)
  config <- yaml::read_yaml(tmpfile)
  
  modelpath <- paste0(home_dir, "montagu/202110gavi-3")
  datapath <- paste0(home_dir, "input_data")
  country <- "AGO"
  nsamples <- 5
  redraw <- TRUE
  random_seed <- 1
  raster_tmp <- ocvImpact::create_incid_raster(modelpath, datapath, country, nsamples, redraw, random_seed)

  ## the file exists
  incid_out_fn <- paste0(datapath, "/incidence/", country, "_incid_5k_", nsamples, ".tif")
  expect_true(file.exists(incid_out_fn))
  ## it has the correct number of layers 
  expect_equal(raster::nlayers(raster_tmp), nsamples)
  ## it has non-NA values 
  expect_true(!all(is.na(raster::values(raster_tmp))))
  ## delete the temporary raster 
  file.remove(incid_out_fn)


  # BGD doesn't work because we didn't upload the raster file for it
  country <- "BGD"
  expect_error(ocvImpact::create_incid_raster(modelpath, datapath, country, nsamples, redraw, random_seed))


  # India works because we use the BGD's IR for India 
  country <- "IND"
  raster_tmp <- ocvImpact::create_incid_raster(modelpath, datapath, country, nsamples, redraw, random_seed)

  ## the file exists
  incid_out_fn <- paste0(datapath, "/incidence/", country, "_incid_5k_", nsamples, ".tif")
  expect_true(file.exists(incid_out_fn))
  ## it has the correct number of layers 
  expect_equal(raster::nlayers(raster_tmp), nsamples)
  ## it has non-NA values 
  expect_true(!all(is.na(raster::values(raster_tmp))))
  ## each layer will only have the same value for all the grid cells 
  expect_equal(dim(table(raster::values(raster_tmp))), nsamples)
  ## delete the temporary raster 
  file.remove(incid_out_fn)

})
