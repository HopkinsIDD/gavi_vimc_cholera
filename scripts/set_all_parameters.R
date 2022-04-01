runname <- "202110gavi-3" ### Most important
scenarios <- c("campaign-default", "no-vaccination")
# countries <- c("COD", "ETH", "KEN", "SOM", "SSD", "HTI")
countries <-    c("AGO", "BDI", "BEN", "BFA", "CAF", "CIV", "CMR", "COD", "COG", "DZA", "ETH", "GHA", "GIN", "GNB", "KEN", "LBR", "MDG", "MLI",
                  "MOZ", "MRT", "MWI", "NAM", "NER", "NGA", "RWA", "SEN", "SLE", "SOM", "SSD", "TCD", "TGO", "TZA", "UGA", "ZAF", "ZMB", "ZWE",
                  "AFG", "HTI", "IRN", "IRQ", "NPL", "PAK", "PHL", "THA", "YEM", "IND", "BGD")
# countries <- c("COD", "ETH", "BGD", "HTI", "YEM")
clean_outputs <- TRUE #changed on 1.31
targeting_strategy <- "affected_pop"
clean_incid <- TRUE #changed on 1.31
num_samples <- 50
num_skip_years <- 3

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
# newly added on 11/18/2021 -- like this for now, but will change to the following
# c("CMR", "COD", "ETH", "IRQ", "KEN", "LBR", "NGA", "SEN", "SOM", "YEM", "ZAF", "ZWE") -- FALSE
# c("HTI", "BGD", "IND") -- not sure


no_country_list <- c("CMR", "COD", "ETH", "IRQ", "KEN", "LBR", "NGA", "SEN", "SOM", "YEM", "ZAF", "ZWE")
not_sure_country_list <- c("HTI", "BGD", "IND") 
use_country_incid_trend <- lapply(countries, function(country_code){
  country_index <- match(country_code, countries)
  use_country_incid_trend[country_index] <- !(country_code %in% no_country_list | country_code %in% not_sure_country_list)
})
use_country_incid_trend <- unlist(use_country_incid_trend)

# rm(country_code)
# rm(country_index)
rm(no_country_list)
rm(not_sure_country_list)

incidence_rate_trend <- TRUE
outbreak_multiplier <- TRUE
### For now, the random seed equals setting number
random_seed <- dplyr::case_when(
  !incidence_rate_trend & !outbreak_multiplier ~ 1, 
  !incidence_rate_trend & outbreak_multiplier ~ 2, 
  incidence_rate_trend & !outbreak_multiplier ~ 3, 
  incidence_rate_trend & outbreak_multiplier ~ 4
)