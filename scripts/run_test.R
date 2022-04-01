#################################
# Author: Kaiyue Zou            #
# Date: 6/29/2021               #
# Purpose: Test Run             #
#################################

### montagu server handle
setwd('C:/Users/ZOU/Desktop')
source("montagu_handle.R")
montagu::montagu_diseases()

### Setup

##new rtools package
#writeLines('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron")
#Sys.which("make")
#install.packages("jsonlite", type = "source")
#install.packages('rtools40', type = "source")

##load the ocvImpact package
setwd('C:/Users/ZOU/Desktop/gavi_vimc_cholera')
roxygen2::roxygenise("packages/ocvImpact")
install.packages("packages/ocvImpact", type = "source", repos = NULL)

### Try to use the function in the package (coverage data)
modelpath = 'C:/Users/ZOU/Desktop/gavi_vimc_cholera/montagu/201910gavi-5' #does extra / make a difference? No. 
ocvImpact::retrieve_montagu_coverage(modelpath)
ocvImpact::retrieve_montagu_coverage(modelpath, scenarios_all = FALSE)

### Try retrieve_montagu_centralburden_template
ocvImpact::retrieve_montagu_centralburden_template(modelpath)
ocvImpact::retrieve_montagu_centralburden_template(modelpath, expectations_all = FALSE)

### Try retrieve_montagu_population
ocvImpact::retrieve_montagu_population(modelpath)
ocvImpact::retrieve_montagu_population(modelpath, countries_all = FALSE)
ocvImpact::retrieve_montagu_population(modelpath, expectations_all = FALSE)
ocvImpact::retrieve_montagu_population(modelpath, expectations_all = FALSE, countries_all = FALSE)

### Try retrieve_montagu_agePop
ocvImpact::retrieve_montagu_agePop(modelpath)

### Try retrieve_montagu_lifeExpectancy
ocvImpact::retrieve_montagu_lifeExpectancy(modelpath)

### Try the import function
country = 'COD'
scenario = 'cholera-campaign-default-test'
ocvImpact::import_coverage_scenario(modelpath, country, scenario)

### Try to get population data for a country
ReturnedPop = ocvImpact::import_country_population(modelpath, country)

### Before run
datapath = "C:/Users/ZOU/Desktop/gavi_vimc_cholera"
nsamples = 30
modelpath = 'C:/Users/ZOU/Desktop/gavi_vimc_cholera/montagu/201910gavi-5'
country = "AFG"
mpathname = modelpath
dpathname = datapath


#*******************************************************************************
### Before Mac Run
setwd('/Users/kaiyuezou/Documents/JHU_Projects/vimc_model/gavi_vimc_cholera')
modelpath <- '/Users/kaiyuezou/Documents/JHU_Projects/vimc_model/gavi_vimc_cholera/montagu/201910gavi-5'
datapath <- '/Users/kaiyuezou/Documents/JHU_Projects/vimc_model/gavi_vimc_cholera/input_data'
nsamples = 30
country = "COD"
mpathname = modelpath
dpathname = datapath

#======Initialize Montagu package======#
if (!require('montagu', character.only = T)) {
  drat:::add("vimc")
  install.packages('montagu')
  library('montagu', character.only = T)
}
source("scripts/montagu_handle.R")

#======Use the ocvImpact package======#
roxygen2::roxygenise("packages/ocvImpact")
install.packages("packages/ocvImpact", type = "source", repos = NULL)
library('ocvImpact', character.only = T)

### Try the new function
library(dplyr)
incid_trend_function <- ocvImpact::generate_flatline_multiplier(trendtype = 'incidence rate', 
                                                                datapath = datapath, 
                                                                modelpath = modelpath, 
                                                                country = country)
incid_trend_function(year = 2050)
year <- 2000:2100
for (years in year){
  if(exists('multiplier')){
    multiplier <- c(multiplier, incid_trend_function(year = years))
  }else{
    multiplier <- incid_trend_function(year = years)
  }
}
plot(year, multiplier, 
     xlab = 'Years', ylab = 'multiplier (2014 as ref)', 
     title(main = country, sub = 'Whether the country own data is used'))


### The for loop
countries <- montagu::montagu_expectation_countries("JHU-Lee", "202110gavi-2", 378)
year <- 2000:2100
datathreshold = 5

pdf('/Users/kaiyuezou/Desktop/CountryIncidenceRate.pdf', width = 5, height = 5)
for (country in countries$id) {
  
  if (country %in% casedf_pop$country_code 
    & length(unique(casedf_pop[casedf_pop$country_code == country, ]$year)) >= as.numeric(datathreshold)) {
    subtitle = 'This country uses its own data in the WHO dataset'
  }else{
    subtitle = 'This country uses the average data from other countries'
  }
  
  incid_trend_function <- ocvImpact::generate_flatline_multiplier(trendtype = 'incidence rate', 
                                                                  datapath = datapath, 
                                                                  modelpath = modelpath, 
                                                                  country = country)
  rm(multiplier)
  for (years in year){
    if(exists('multiplier')){
      multiplier <- c(multiplier, incid_trend_function(year = years))
    }else{
      multiplier <- incid_trend_function(year = years)
    }
  }
  plot(year, multiplier, 
       xlab = 'Years', ylab = 'multiplier (year 2014 as ref)', 
       title(main = country, sub = subtitle))
  
}
dev.off()


#*******************************************************************************
### Try out the different outbreak multiplier
outbreak_2000 <- outbreak_trend_function(1)
raster::plot(raster::raster(outbreak_2000, layer = 1))
outbreak_2014 <- outbreak_trend_function(15)
raster::plot(raster::raster(outbreak_2014, layer = 5))
outbreak_trend_function(101)

pdf('/Users/kaiyuezou/Desktop/Outbreak_Multiplier_plotting3.pdf')
for (i in 1:101) {
  year <- 1999 + i
  outbreak_one_year <- outbreak_trend_function(i)
  raster::plot(raster::raster(outbreak_one_year, layer = 1), 
               main = paste0('Year ', year, ' outbreak multiplier for COD'), 
               sub = 'The first layer out of 30 layers')
}
dev.off()

