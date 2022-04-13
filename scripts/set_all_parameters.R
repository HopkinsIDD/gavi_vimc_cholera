#====== The Basics ======#
runname <- "202110gavi-3" ### Most important -- this is borrowed for the new surveillance project to pull demo data easily from Montagu
scenarios <- c("campaign-default", "no-vaccination")
countries <-c("AGO", "BDI", "BEN", "BFA", "CAF", "CIV", "CMR", "COD", "COG", "DZA", "ETH", "GHA", "GIN", "GNB", "KEN", "LBR", "MDG", "MLI",
              "MOZ", "MRT", "MWI", "NAM", "NER", "NGA", "RWA", "SEN", "SLE", "SOM", "SSD", "TCD", "TGO", "TZA", "UGA", "ZAF", "ZMB", "ZWE",
              "AFG", "HTI", "IRN", "IRQ", "NPL", "PAK", "PHL", "THA", "YEM", "IND", "BGD")
countries <- c("COD", "ETH", "BGD", "HTI", "YEM") #for testing, temporary 

#====== Surveillance Project Specific ======#
targeting_strategy <- "surveillance_project" #c("surveillance_project", "affected_pop", "incidence")
vac_incid_threshold <- 1/1000 #c(1/1000, 1/5000, 1/10000), use one at a time for now
vac_unconstrained <- TRUE
vac_admin_level <- 2 #c(1, 2)
vac_coverage <- 0.8 #for now
surveillance_scenarios <- c("no-estimate", "global-estimate", "district-estimate") #use all three at the same time for now

#====== Other Parameters ======#
num_skip_years <- 3 #relevant to the surveillance project
num_samples <- 50 #relevant to the surveillance project
clean_outputs <- TRUE #changed on 1.31
clean_incid <- TRUE #changed on 1.31

# ### tmp
# outbreak_file_name <- paste0("/home/kaiyuezou/VIMC_Model/202110gavi-3/gavi_vimc_cholera/input_data", '/outbreak/outbreak_df.csv')
# if (file.exists(outbreak_file_name)){
#   outbreak_df <- readr::read_csv(outbreak_file_name)
# }else{
#   stop('The multi-country outbreak dataset does not exist under the correct directory. Please check. ')
# }
# outbreak_countries <- unique(outbreak_df$country)[!unique(outbreak_df$country) %in% c("TZA_zanzibar")]
# # countries <- countries[!countries %in% outbreak_countries]
# countries <- outbreak_countries

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

#====== Setting Parameters--incidence rate trend and outbreak multiplier ======#
incidence_rate_trend <- FALSE #for surveillance for now
outbreak_multiplier <- FALSE #for surveillance for now
### For now, the random seed equals setting number
random_seed <- dplyr::case_when(
  !incidence_rate_trend & !outbreak_multiplier ~ 1, 
  !incidence_rate_trend & outbreak_multiplier ~ 2, 
  incidence_rate_trend & !outbreak_multiplier ~ 3, 
  incidence_rate_trend & outbreak_multiplier ~ 4
)
