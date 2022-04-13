# Function to update incidence column based on the ec_list (expected cases rasterStack)
update_incidence <- function(rc_list, # made into sf so that can extract based on shapefiles
                             ec_list
                             ){
  
  next_year <- max(rc_list[[1]]$target_year)
  message(paste0("Updating incidence in preparation of next year's (", next_year, ") campaign."))

  rc_this_year_admin1 <- rc_list[[1]] %>% filter(target_year == next_year)
  rc_this_year_admin2 <- rc_list[[2]] %>% filter(target_year == next_year)
  
  yr_index <- which(next_year-1 == unique(rc_list[[1]]$target_year))
  ec_this_year_admin1 <- raster::subset(ec_list[["ec_rasterStack_admin1"]], yr_index, drop = FALSE)
  ec_this_year_admin2 <- raster::subset(ec_list[["ec_rasterStack_admin2"]], yr_index, drop = FALSE)
  
  # get the expected case number of each places at admin1 and admin2 level
  ec_admin1 <- exactextractr::exact_extract(ec_this_year_admin1, rc_this_year_admin1, "sum")
  ec_admin2 <- exactextractr::exact_extract(ec_this_year_admin2, rc_this_year_admin2, "sum")
  
  # get the incidence
  new_inc_admin1 <- ec_admin1 / rc_this_year_admin1$pop_model
  new_inc_admin2 <- ec_admin2 / rc_this_year_admin2$pop_model
  
  # update tables in rc_list
  rc_this_year_admin1$incidence <- new_inc_admin1
  rc_this_year_admin2$incidence <- new_inc_admin2
  
  rc_list[[1]] <- rc_list[[1]] %>% filter(target_year < next_year) %>% rbind(rc_this_year_admin1)
  rc_list[[2]] <- rc_list[[2]] %>% filter(target_year < next_year) %>% rbind(rc_this_year_admin2)
  
  rm(ec_this_year_admin1, ec_this_year_admin2)
  gc()
  
  return(rc_list)
  
}
