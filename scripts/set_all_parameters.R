#====== The Basics ======#
targeting_strategy <- "affected_pop" #c("threshold_unconstrained", "affected_pop", "incidence"), "threshold_unconstrained" means it's for the surveillance project
#runname <- ifelse(targeting_strategy == "threshold_unconstrained", "202302_survms", "202110gavi-3") ### Most important -- this is borrowed for the new surveillance project to pull demo data easily from Montagu
runname <- ifelse(targeting_strategy == "threshold_unconstrained", "202302_survms", "202310gavi-4") ### for 2023 VIMC core:Most important -- this is borrowed for the new surveillance project to pull demo data easily from Montagu


#====== Shared parameters ======#
scenarios <- c("campaign-default", "no-vaccination")
num_skip_years <- 3   #district-level skipped years, relevant to both projects
num_samples <- 100    #shared by both projects
use_random_seed <- TRUE   #whether or not to have a random seed that governs the stochasticity

if (runname == "202310gavi-4"){
  self_random_seed <- 33  #use the same random seed for all setting and scenarios for the 202310gavi-4 touchstone
} else {
  self_random_seed <- NULL  #for now, just use the random seed specified by the setting
}

clean_outputs <- TRUE
clean_incid <- FALSE
if (runname == "202310gavi-4"){
  ndoses <- c("one", "two") ##calam added option for the number of doses, only applies to the VIMC Core model 2023 touchstone
}
# default country list
ids <- readr::read_csv("input_data/locations_todeletelater.csv") # location ids
cw <- readr::read_csv("input_data/region_country.csv")
locs <- dplyr::right_join(cw, ids, by = c("country" = "region")) ## region & country grouping
locs <- dplyr::filter(locs, !is.na(region))
countries <- unique(locs$country)
# countries simulated in the VIMC core project
countries <-c("AGO", "BDI", "BEN", "BFA", "CAF", "CIV", "CMR", "COD", "COG", "DZA", "ETH", "GHA", "GIN", "GNB", "KEN", "LBR", "MDG", "MLI",
              "MOZ", "MRT", "MWI", "NAM", "NER", "NGA", "RWA", "SEN", "SLE", "SOM", "SSD", "TCD", "TGO", "TZA", "UGA", "ZAF", "ZMB", "ZWE",
              "AFG", "HTI", "IRN", "IRQ", "NPL", "PAK", "PHL", "THA", "YEM", "IND", "BGD") #will likely to only include the countries in sub-Saharan Africa
# countries simulated in the surveillance project              
countries <-c("AGO", "BDI", "BEN", "BFA", "CAF", "CIV", "CMR", "COD", "COG", "DZA", "ETH", "GHA", "GIN", "GNB", "KEN", "LBR", "MDG", "MLI",
              "MOZ", "MRT", "MWI", "NAM", "NER", "NGA", "RWA", "SEN", "SLE", "SOM", "SSD", "TCD", "TGO", "TZA", "UGA", "ZAF", "ZMB", "ZWE") # now only includes the countries in Africa


#====== Surveillance Project Specific ======#
save_intermediate_raster <- TRUE #if TRUE, this will save memory during the model run
save_final_output_raster <- TRUE #if TRUE, all the raster files will be saved
ir_pre_screening <- FALSE #if TRUE, countries with low ir where no vacc will happen will be skipped 
vac_incid_thresholds <- c(1/1000, 1/5000, 1/10000) 
vac_unconstrained <- TRUE #or refer to an external coverage dataset 
vac_admin_level <- "both" #c("both", "admin1", "admin2"), running both admin levels for now, running a single one when issues occur 
vac_coverage <- 0.68 #borrowed from the method section of the previous study (https://doi.org/10.1371/journal.pmed.1003003) 
surveillance_scenarios <- c("no-estimate", "global-estimate", "district-estimate") #use all three at the same time for now; the estimates can be drawn from the external dataset instead of being teemed in the config 
testing_sensitivity <- list("no-estimate" = 1, "global-estimate" = 0.8, "district-estimate" = 1)
vac_interval <- 1 #country level, unit: years
sim_start_year <- 2022
vac_start_year <- 2022 #we assume that vaccination campaign starts and ends at the beginning of each year 
vac_end_year <- 2030
sim_end_year <- 2035 
use_mean_ir <- TRUE
mean_ir_span <- 5


#====== Shared parameters -- related to country-level incidence rate trends ======#
use_country_incid_trend <- rep(TRUE, length(countries))

no_country_list <- c("CMR", "COD", "ETH", "IRQ", "KEN", "LBR", "NGA", "SEN", "SOM", "YEM", "ZAF", "ZWE")
not_sure_country_list <- c("HTI", "BGD", "IND") 
use_country_incid_trend <- lapply(countries, function(country_code){
  country_index <- match(country_code, countries)
  use_country_incid_trend[country_index] <- !(country_code %in% no_country_list | country_code %in% not_sure_country_list)
})
use_country_incid_trend <- unlist(use_country_incid_trend)

rm(no_country_list)
rm(not_sure_country_list)


#====== Setting Parameters -- based on incidence rate trend and outbreak multiplier ======#
incidence_rate_trend <- TRUE
outbreak_multiplier <- FALSE 
# For now, the random seed equals setting number
if(use_random_seed & is.null(self_random_seed)){
  random_seed <- dplyr::case_when(
    !incidence_rate_trend & !outbreak_multiplier ~ 1, 
    !incidence_rate_trend & outbreak_multiplier ~ 2, 
    incidence_rate_trend & !outbreak_multiplier ~ 3, 
    incidence_rate_trend & outbreak_multiplier ~ 4
  )
}else if(use_random_seed & !is.null(self_random_seed)){
  random_seed <- self_random_seed
}else if(!use_random_seed){
  random_seed <- sample(1:10000, 1) #let loose on random seed 
}
