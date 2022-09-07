#' @name check_table_screening
#' @title check_table_screening
#' @description check_table_screening
#' @param spathname path to input incidence data 
#' @param country country code 
#' @param scenario
#' @param threshold incidence rate threshold
#' @param vac_admin_level the admin level to check
#' @return 
#' @export
#' @include
check_table_screening <- function(spathname, country, scenario, threshold, vac_admin_level){
  
  scr_fn <- paste(spathname, "country_simulation_skipped.csv", sep = "/")
    
  if(!file.exists(scr_fn)){
    return(FALSE)
  
  }else{
    scr_df <- readr::read_csv(scr_fn)
    if(country %in% scr_df$country_skipped){
      if(scr_df[scr_df$country_skipped == country, ]$knowable){
        if(vac_admin_level == "both" & (any(scr_df[scr_df$country_skipped == country, 
                                                  c("admin1_highest_incidence_rate", "admin2_highest_incidence_rate")]) >= as.numeric(threshold))){
          cache$novacc_campde_transfer <- FALSE
          return(TRUE)
        }else if(vac_admin_level == "admin1" & (scr_df[scr_df$country_skipped == country, c("admin1_highest_incidence_rate")] >= as.numeric(threshold))){
          cache$novacc_campde_transfer <- FALSE
          return(TRUE)
        }else if(vac_admin_level == "admin2" & (scr_df[scr_df$country_skipped == country, c("admin2_highest_incidence_rate")] >= as.numeric(threshold))){
          cache$novacc_campde_transfer <- FALSE
          return(TRUE)
        }else{
          if(scenario == "campaign-default"){
            message(paste(" -- Country", country, "will be skipped fully for the current simulation scenario due to its low incidence rate. "))
            novacc_campde_transfer(cache$rawoutpath, cache$config)
            quit()
          }else if(scenario == "no-vaccination"){
            message(paste(" -- Country", country, "will only be simulated for the current no-vaccination scenario due to its low incidence rate. "))
            cache$novacc_campde_transfer <- TRUE
            return(TRUE)
          }
          
        }
      }else{return(TRUE)}
    
    }else{return(FALSE)}
  }
}

#' @name update_table_screening
#' @title update_table_screening
#' @description update_table_screening
#' @param datapath path to the input folder
#' @param modelpath path to montagu
#' @param country country code 
#' @param scenario 
#' @param threshold current incidence rate threshold
#' @param vac_start_year the year when the vacc is started
#' @param vac_end_year the year when the vacc ends 
#' @param rc_list a list of targeting table to check confirmed_incidence_rate
#' @param rc_targeted admin levels to be targeted at 
#' @param nsamples number of layers 
#' @return 
#' @export
#' @include
update_table_screening <- function(datapath, modelpath, country, scenario, threshold, vac_start_year, vac_end_year, rc_list, rc_targeted, nsamples){
  # Get the dataset
  scr_fn <- paste(datapath, "incidence", "country_simulation_skipped.csv", sep = "/")
  if(!file.exists(scr_fn)){
    scr_df <- tibble::tibble( country_skipped = as.character(), admin1_highest_incidence_rate = as.numeric(), 
                              admin2_highest_incidence_rate = as.numeric(), knowable = as.logical())
  }else{
    scr_df <- readr::read_csv(scr_fn)
  }
  
  # Check the population change
  pop_df <- ocvImpact::import_country_population(modelpath, country, redownload = FALSE)
  pop_prior <- 0
  for (years in vac_start_year:vac_end_year){
    pop <- pop_df[pop_df$year == years, ]$pop_model
    if(pop <- pop_prior){
      scr_df <- scr_df %>% add_row(country_skipped = country, admin1_highest_incidence_rate = NA, admin2_highest_incidence_rate = NA, knowable = FALSE)
      readr::write_csv(scr_df, scr_fn)
      message(paste(" -- Country", country, "will be targeted at least once for the current simulation scenario, full simulation will begin. "))
      return()
    }
    pop_prior <- pop
  }
  
  # Get the max of the confirmed_incidence_rate
  cir1 <- c()
  cir2 <- c()
  for(layer_idx in 1:nsamples){
    if("rc1" %in% rc_targeted){cir1 <- c( cir1, max(rc_list[[layer_idx]]$rc1$confirmed_incidence_rate) )}
    if("rc2" %in% rc_targeted){cir2 <- c( cir2, max(rc_list[[layer_idx]]$rc2$confirmed_incidence_rate) )}
  }
  
  # Record
  if(!country %in% scr_df$country_skipped){
    scr_df <- scr_df %>% 
      add_row(country_skipped = country, admin1_highest_incidence_rate = ifelse(is.infinite(max(cir1)), NA, max(cir1)), 
              admin2_highest_incidence_rate = ifelse(is.infinite(max(cir2)), NA, max(cir2)), knowable = TRUE)
  }else{
    if("rc1" %in% rc_targeted){scr_df[scr_df$country_skipped == country, ]$admin1_highest_incidence_rate <- max(cir1)}
    if("rc2" %in% rc_targeted){scr_df[scr_df$country_skipped == country, ]$admin2_highest_incidence_rate <- max(cir2)}
  }
  readr::write_csv(scr_df, scr_fn)
  
  # Compare
  if(all(c("rc1", "rc2") %in% rc_targeted)){ vac_admin_level <- "both"
    }else{ vac_admin_level <- c("admin1", "admin2")[grepl(stringr::str_extract(rc_targeted, "[0-9]{1}"), c("admin1", "admin2"))] }       
  cache$ir_pre_screening_pass <- check_table_screening(spathname = file.path(datapath, "incidence"), country, scenario, threshold, vac_admin_level)

}

#' @name novacc_campde_transfer 
#' @title novacc_campde_transfer 
#' @description novacc_campde_transfer 
#' @param rawoutpath path to the raw output folder 
#' @param config the config with all the parameters 
#' @return 
#' @export
#' @include
novacc_campde_transfer <- function(rawoutpath, config){
  
  message(" *----* Using the novacc scenario model output for the camp scenario. ")
  incidence_rate_trend <- config$setting$incidence_rate_trend
  outbreak_multiplier <- config$setting$outbreak_multiplier
  vac_incid_threshold <- config$vac$vac_incid_threshold
  surveillance_scenarios <- c("no-estimate", "global-estimate", "district-estimate")
  country <- config$country

  for(surveillance_scenario in surveillance_scenarios){
    rc1_out_fn_novacc <- paste0(rawoutpath, "/no-vaccination/", paste("incid", incidence_rate_trend, "outbk", outbreak_multiplier, 
                          vac_incid_threshold, surveillance_scenario, country, sep = "_"), "_rc_admin1.csv")
    rc2_out_fn_novacc <- paste0(rawoutpath, "/no-vaccination/", paste("incid", incidence_rate_trend, "outbk", outbreak_multiplier, 
                          vac_incid_threshold, surveillance_scenario, country, sep = "_"), "_rc_admin2.csv")
    rc1_out_fn_campde <- paste0(rawoutpath, "/campaign-default/", paste("incid", incidence_rate_trend, "outbk", outbreak_multiplier, 
                          vac_incid_threshold, surveillance_scenarios, country, sep = "_"), "_rc_admin1.csv")
    rc2_out_fn_campde <- paste0(rawoutpath, "/campaign-default/", paste("incid", incidence_rate_trend, "outbk", outbreak_multiplier, 
                          vac_incid_threshold, surveillance_scenarios, country, sep = "_"), "_rc_admin2.csv")
    if(file.exists(rc1_out_fn_novacc)){
      if(any(file.exists(rc1_out_fn_campde))){
        ow_idc <- !length(unique(tools::md5sum(c(rc1_out_fn_novacc, rc1_out_fn_campde[file.exists(rc1_out_fn_campde)])))) == 1 
      }else{ow_idc <- TRUE}
      file.copy(from = rc1_out_fn_novacc, to = rc1_out_fn_campde, overwrite = ow_idc)
    }  
    if(file.exists(rc2_out_fn_novacc)){
      if(any(file.exists(rc2_out_fn_campde))){
        ow_idc <- !length(unique(tools::md5sum(c(rc2_out_fn_novacc, rc2_out_fn_campde[file.exists(rc2_out_fn_campde)])))) == 1 
      }else{ow_idc <- TRUE}
      file.copy(from = rc2_out_fn_novacc, to = rc2_out_fn_campde, overwrite = ow_idc)
    } 
  }
  
}
