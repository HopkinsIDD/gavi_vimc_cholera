#################################
# Author: Kaiyue Zou            #
# Date: 6/29/2021               #
# Purpose: Test Run             #
#################################

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
modelpath = 'C:/Users/ZOU/Desktop/only local ocvImpact/montagu/201910gavi-5' #does extra / make a difference? No. 
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







