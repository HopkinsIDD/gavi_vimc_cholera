## Script to process the "custom" health zone vector for the DRC Case study to match the gadm admin 2 vector 

##libraries
library(sf)
library(readr)
library(tidyr)
library(dplyr)

##load the gadm vector used in the VIMC Core Model
gadm_vector <- readRDS("COD_adm2.sf.rds")
#class(gadm_vector)
#plot(gadm_vector[1])

##load health zone vector from humanitarian data exchange
custom_vector <- sf::st_read("osm_rdc_sante_zones_211212.gpkg")
#class(custom_vector)
#plot(custom_vector[1])

##Procedure 1: check names from the DRC Targeting table match those inside the DRC Health Zone vector 
## This procedure is optional, use to inspect the DRC Health zone vector

##load DRC health zone names from the DRC targeting table


##DRC_targeting_raw <- readr::read_csv("DRC_targeting_raw.csv", 
                              ##col_types = cols(target_population = col_skip(), 
                                               ##Vaccinated_people_2_rounds = col_skip(), 
                                               ##gadm_match = col_skip(), Notes = col_skip()))

##check_names <- which(DRC_targeting_raw$Health_zone %in% custom_vector$name)
##match_names <- DRC_targeting_raw[check_names,] ## targeting table health zones with an exact name match to the custom vector
##no_match_names <- DRC_targeting_raw[-check_names,] ## targeting table health zones with no match


## Procedure 2: subset the custom vector to keep only health zones that match the targeting table
## This procedure is optional, use to get a vector with the DRC Health zones mentioned in the DRC targeting table

##check_names_vector <- which(custom_vector$name %in% DRC_targeting_raw$Health_zone)
##health_zones_needed <- custom_vector[check_names_vector,]
##plot(health_zones_needed[1])


## Procedure 3
## This procedure is necessary in order to generate a health zone vector that is compatible with our pipeline

##make custom_vector match names, resolution, crs, object class of gadm shapefile

##rename and add columns
custom_vector_final <- custom_vector %>%
  dplyr::rename(NAME_2 = name, geometry = geom) %>%
  dplyr::select(NAME_2, geometry) %>%
  dplyr::mutate(GID_0 = "COD", NAME_0 = "Democratic Republic of the Congo", GID_2 = "<NA>", NL_NAME_1 = "<NA>", VARNAME_2 = "<NA>", NL_NAME_2 = "<NA>", TYPE_2 = "<NA>", 
                ENGTYPE_2 = "<NA>", CC_2 = "<NA>", HASC_2  = "<NA>")


## get GID_1 and NAME_1 (province data) from the gadm admin1 shapefile with a spatial join


## get admin 1 gadm data

gadm_admin1 <- sf::st_as_sf(geodata::gadm("COD", level = 1, path = "."))

## get GID_1 and NAME_1 (province data) from the gadm shapefile with a spatial join

sf1_sf2_join <- sf::st_join(custom_vector_final, gadm_admin1, join = st_nearest_feature)

## extract the values of "GID_1" from the gadm vector for locations in the custom vector
custom_vector_final$GID_1 <- sf1_sf2_join$GID_1

## get NAME_1 values
custom_vector_final$NAME_1 <- sf1_sf2_join$NAME_1

## remove joined to save memory
rm(sf1_sf2_join)

## rename health zones to match our DRC targeting table

custom_vector_final <- custom_vector_final %>%
  dplyr:: mutate(NAME_2 = dplyr::case_when(
    NAME_2 == "Nyirangongo" ~ "Nyiragongo",
    NAME_2 == "Miti-Murhesa" ~ "Miti Murhesa",
    NAME_2 == "Kabondo-Dianda" ~ "Kabondo Dianda",
    NAME_2 == "Lolanga Mampoko" ~ "Mampoko",
    TRUE ~ NAME_2  # Keep the original value if none of the conditions are met
  ))

## make sure gadm and custom vectors have the same crs

if (sf::st_crs(gadm_vector) == sf::st_crs(custom_vector_final)){
  message("The gadm admin 2 layer and the custom health zone layer have the same crs")
} else {
  sf::st_crs(custom_vector_final) <- sf::st_crs(gadm_vector)
  message("The custom health zone layer crs has been set to that of the gadm admin 2 layer")
}

## redorder columns according to the column order in the gadm vector

if (all(colnames(gadm_vector) %in% colnames(custom_vector_final))){
  custom_vector_final <- custom_vector_final %>%
    dplyr::select(colnames(gadm_vector))
} else {
  stop("The custom vector and the gadm admin 2 vector do not have the same column names")
}

## write custom vector to file

saveRDS(custom_vector_final, "custom_shapefile.rds")
