### quick raw output check for the surveillance project -- just for .csv files 
dir <- "/..." #this is the directory of the raw output folder 
iso <- stringr::str_extract(list.files(dir), "_[A-Z]{3}_.*.csv")
iso <- iso[!is.na(iso)]
iso <- stringr::str_extract(iso, "[A-Z]{3}")
iso <- sort(table(iso))
iso_sml <- iso[iso == 3]
length(iso_sml)

### quick raw output check for the surveillance project -- all raw output
dir <- "/..." #this is the directory of the raw output folder 
ec <- list.files(dir)[stringr::str_detect(list.files(dir), "_[A-Z]{3}_ec_")]
iso <- stringr::str_extract(ec, "_[A-Z]{3}_")
iso <- stringr::str_extract(iso, "[A-Z]{3}")
bad_ec <- names(iso)[iso != 2*2*14]
print(paste0("Following countries have incomplete expected case raster files: ", bad_ec))

sus <- list.files(dir)[stringr::str_detect(list.files(dir), "_[A-Z]{3}_sus_")]
iso <- stringr::str_extract(sus, "_[A-Z]{3}_")
iso <- stringr::str_extract(iso, "[A-Z]{3}")
bad_sus <- names(iso)[iso != 2*2*14]
print(paste0("Following countries have incomplete susceptibility raster files: ", bad_sus))

vac <- list.files(dir)[stringr::str_detect(list.files(dir), "_[A-Z]{3}_vac_")]
iso <- stringr::str_extract(vac, "_[A-Z]{3}_")
iso <- stringr::str_extract(iso, "[A-Z]{3}")
bad_vac <- names(iso)[iso != 2*2*14]
print(paste0("Following countries have incomplete vaccination raster files: ", bad_vac))
