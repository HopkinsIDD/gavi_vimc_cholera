# roxygen2::roxygenise("packages/ocvImpact")
# install.packages("packages/ocvImpact", type = "source", repos = NULL)

#######Kaiyue Added on 8/15/2021#######
#======Use other packages needed======#
chooseCRANmirror(ind = 77) #specify the mirror so that the packages can be successfully installed in the non-interactive way

package_list <- c(
  "GADMTools", 
  "rgdal", 
  "drat", 
  "roxygen2", 
  "data.table",
  "dplyr",
  "exactextractr",
  "fasterize",
  "optparse",
  "purrr",
  "raster",
  "readr",
  "sf",
  "stringr",
  "readxl",
  "tibble",
  "tidyr",
  "yaml"
)

for (package in package_list) {
  if (!require(package = package, character.only = T)) {
    install.packages(pkgs = package)
    library(package = package, character.only = T)
  }
  detach(pos = which(grepl(package, search())))
}

#======Initialize Montagu package======#
if (!require('montagu', character.only = T)) {
  drat:::add("vimc")
  install.packages('montagu')
  library('montagu', character.only = T)
}
source("scripts/montagu_handle.R")

#======Use the ocvImpact package======#
if (!require('ocvImpact', character.only = T)) {
  roxygen2::roxygenise("packages/ocvImpact")
  install.packages("packages/ocvImpact", type = "source", repos = NULL)
  library('ocvImpact', character.only = T)
}

#======For the convenience of debugging======#
###These a few lines can be deleted safely after the model can run smoothly on the server. 
library(raster)
roxygen2::roxygenise("packages/ocvImpact")
install.packages("packages/ocvImpact", type = "source", repos = NULL)
library('ocvImpact', character.only = T)

###########Comment completed###########

library(magrittr)
source("scripts/set_all_parameters.R") ## read in parameters
opathname <- file.path("output_final", runname)
mpathname <- file.path("montagu", runname)

for (i in 1:length(scenarios)){

  scn <- scenarios[i]

  ## Stochastic outputs
  incidence_rate_trend <-TRUE #just for now
  outbreak_multiplier <- TRUE #just for now
  single_setting <- FALSE
  if(single_setting){
    suffix <- paste0(".*", "_incid_trend_", incidence_rate_trend, "_outb_layer_", outbreak_multiplier, "_stoch.csv")
  }else{
    suffix <- ".*_stoch.csv"
  }
  ##added to specify correct filenames for the outputs of the 202310gavi-4 touchstone
  if (runname == "202310gavi-4"){
    ndoses <- config$vacc$ndoses
    stoch_fns <- list.files(opathname, pattern = paste0(scn,ndoses,suffix))
  } else{
    stoch_fns <- list.files(opathname, pattern = paste0(scn, suffix))
  }
  ##end addition
  stoch_out <- purrr::map_dfr(1:length(stoch_fns), function(j){
    if(single_setting){
      return(readr::read_csv(file.path(opathname, stoch_fns[j])))
    }else{
      stoch <- readr::read_csv(file.path(opathname, stoch_fns[j]))
      filename <- stoch_fns[j]
      ##calam added to find new filenames for the 2023gavi-4 touchstone
      if (runname == "202310gavi-4"){
        inci <- as.logical(stringr::str_split(filename, "_")[[1]][6])
        outb <- as.logical(stringr::str_split(filename, "_")[[1]][9])
      } else {
        inci <- as.logical(stringr::str_split(filename, "_")[[1]][5])
        outb <- as.logical(stringr::str_split(filename, "_")[[1]][8])
      }

      # stoch$run_id <- as.numeric(stoch$run_id)
      stoch <- stoch %>%
        dplyr::mutate(run_id = dplyr::case_when(
          !inci & !outb ~ run_id,
          !inci & outb ~ run_id + 50,
          inci & !outb ~ run_id + 100,
          inci & outb ~ run_id + 150
        ))
      return(stoch)
    }
    
  })
  stoch_final_fn <- paste0(opathname, "/", ocvImpact::import_templateFilename_prefix("stochastic", mpathname), scn, ".csv")
  message(paste("Write stochastic final output:", stoch_final_fn))
  readr::write_csv(stoch_out, stoch_final_fn)

  ## Central outputs
  central_out <- ocvImpact::export_central_template(stoch_out)
  central_final_fn <- paste0(opathname, "/", ocvImpact::import_templateFilename_prefix("central", mpathname), scn, ".csv")
  message(paste("Write central final output:", central_final_fn))
  readr::write_csv(central_out, central_final_fn)

  ## Parameter outputs
  # incidence_rate_trend <- TRUE #just for now
  # outbreak_multiplier <- FALSE #just for now
  # single_setting <- TRUE
  if(single_setting){
    suffix <- paste0(".*", "_incid_trend_", incidence_rate_trend, "_outb_layer_", outbreak_multiplier, "_pars.csv")

    par_fns <- list.files(opathname, pattern = paste0(scn, suffix))
    par_out <- purrr::map_dfr(1:length(par_fns), function(k){
      return(readr::read_csv(file.path(opathname, par_fns[k])))
    }) %>%
      tidyr::pivot_wider(names_from = country, names_sep = ":", values_from = c(aoi, incid_rate, cfr))
    colnames <- stringr::str_split(names(par_out), ":", simplify=TRUE)
    names(par_out) <- stringr::str_remove(paste(colnames[,2], colnames[,1], sep = ":"), pattern = "^:")
    
    par_final_fn <- paste0(opathname, "/", ocvImpact::import_templateFilename_prefix("parameter", mpathname), scn, ".csv")
    message(paste("Write parameter final output:", par_final_fn))
    readr::write_csv(par_out, par_final_fn)

  }else{
    suffix <- c(paste0(".*", "_incid_trend_", "FALSE", "_outb_layer_", "FALSE", "_pars.csv"), 
                paste0(".*", "_incid_trend_", "FALSE", "_outb_layer_", "TRUE", "_pars.csv"), 
                paste0(".*", "_incid_trend_", "TRUE", "_outb_layer_", "FALSE", "_pars.csv"), 
                paste0(".*", "_incid_trend_", "TRUE", "_outb_layer_", "TRUE", "_pars.csv") )

    for(suff in suffix){
      par_fns <- list.files(opathname, pattern = paste0(scn, suff))
      par_out <- purrr::map_dfr(1:length(par_fns), function(k){
        return(readr::read_csv(file.path(opathname, par_fns[k])))
      }) %>%
        tidyr::pivot_wider(names_from = country, names_sep = ":", values_from = c(aoi, incid_rate, cfr))
      colnames <- stringr::str_split(names(par_out), ":", simplify=TRUE)
      names(par_out) <- stringr::str_remove(paste(colnames[,2], colnames[,1], sep = ":"), pattern = "^:")

      if(!exists("scenario_index")){
        scenario_index <- i
      }else{
        if(scenario_index != i){
          rm(par_out_total) #start with a new scenario
          scenario_index <- i
        }

      }

      if(exists("par_out_total")){
        par_out_total <- rbind(par_out_total, par_out)
      }else{
        par_out_total <- par_out
      }
      
    }
    
    par_final_fn <- paste0(opathname, "/", ocvImpact::import_templateFilename_prefix("parameter", mpathname), scn, "-all_setting.csv")
    message(paste("Write parameter final output:", par_final_fn))
    par_out_total[, 1] <- 1:nrow(par_out_total) #run id fix
    readr::write_csv(par_out_total, par_final_fn)
  }

  # par_fns <- list.files(opathname, pattern = paste0(scn, suffix))
  # par_out <- purrr::map_dfr(1:length(par_fns), function(k){
  #   return(readr::read_csv(file.path(opathname, par_fns[k])))
  # }) %>%
  #   tidyr::pivot_wider(names_from = country, names_sep = ":", values_from = c(aoi, incid_rate, cfr))
  # colnames <- stringr::str_split(names(par_out), ":", simplify=TRUE)
  # names(par_out) <- stringr::str_remove(paste(colnames[,2], colnames[,1], sep = ":"), pattern = "^:")
  
  # par_final_fn <- paste0(opathname, "/", ocvImpact::import_templateFilename_prefix("parameter", mpathname), scn, ".csv")
  # message(paste("Write parameter final output:", par_final_fn))
  # readr::write_csv(par_out, par_final_fn)
}
