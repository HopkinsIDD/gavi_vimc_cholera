## aggregate to 20k for alignment

## finally decided we needed 1x1km raster stack to ensure best alignment
pop1k <- raster::stack("input_data/worldpop/ppp_2020_1km_Aggregated.tif")

pop20k <- raster::aggregate(pop1k, 5, filename = "input_data/worldpop/ppp_2020_5km_Aggregated.tif")
message("End script")
