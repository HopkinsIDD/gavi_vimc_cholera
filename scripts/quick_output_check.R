### quick output check
dir <- "/..." #this is the directory of the raw output folder 
iso <- stringr::str_extract(list.files(dir), "_[A-Z]{3}_.*.csv")
iso <- iso[!is.na(iso)]
iso <- stringr::str_extract(iso, "[A-Z]{3}")
iso <- sort(table(iso))
iso_sml <- iso[iso == 3]
length(iso_sml)
