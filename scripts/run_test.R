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
setwd('C:/Users/ZOU/Desktop/ocvImpact new edit')
roxygen2::roxygenise("packages/ocvImpact")
install.packages("packages/ocvImpact", type = "source", repos = NULL)

### Try to use the function in the package (coverage data)
modelpath = 'C:/Users/ZOU/Desktop/ocvImpact new edit/montagu/201910gavi-5/' #does extra / make a difference? No. 
touchstone = '201910gavi-5'
ocvImpact::retrieve_montagu_coverage(modelpath, touchstone)
ocvImpact::retrieve_montagu_coverage(modelpath, touchstone, scenarios_all = FALSE)

### Try retrieve_montagu_centralburden_template
ocvImpact::retrieve_montagu_centralburden_template(modelpath, touchstone)
ocvImpact::retrieve_montagu_centralburden_template(modelpath, touchstone, expectations_all = FALSE)

### Try retrieve_montagu_population
ocvImpact::retrieve_montagu_population(modelpath, touchstone)
ocvImpact::retrieve_montagu_population(modelpath, touchstone, countries_all = FALSE)
ocvImpact::retrieve_montagu_population(modelpath, touchstone, expectations_all = FALSE)
ocvImpact::retrieve_montagu_population(modelpath, touchstone, expectations_all = FALSE, countries_all = FALSE)

### Try retrieve_montagu_agePop
ocvImpact::retrieve_montagu_agePop(modelpath, touchstone)

### Try retrieve_montagu_lifeExpectancy
ocvImpact::retrieve_montagu_lifeExpectancy(modelpath, touchstone)


