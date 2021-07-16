
#' @name get_singular_estimate
#' @title get_singular_estimate
#' @description Generate a dataset with singular estimate for non-raster countries
#' @param datapath path to data 
#' @param country country code
#' @param nyears numeric, number of years we want to calculate average annual #cases from (default = 5 and always most recent available)
#' @return an updated .csv file
#' @export

get_singular_estimate <- function(datapath, country, nyears = 5){
  
  #First decide whether we need to generate singular incidence estimate for this country
  IncidenceTable <- read.csv(file = paste0(datapath, '/incidence/VIMC-47-countries-for-cholera-modelling.csv'))
  NonRasterCountry <- subset(IncidenceTable, !is.na(IncidenceTable$singular_estimate))$country_code
  if (country %in% NonRasterCountry){
    #Read in the WHO data source and calculate
    WHOdf <- subset(read.csv(file = paste0(datapath, '/incidence/who_case_repo_source.csv')), ISO3 == country)
    if (length(WHOdf$year) >= nyears){
      YearList <- list(sort(WHOdf$year, decreasing = TRUE)[1:nyears])
    } else{
      stop(paste("The WHO data source doesn't have ", nyears, " of incidence data for ", country, ". "))
    }
    AverageAnnualCases <- mean(WHOdf$sCh[match(YearList[[1]], WHOdf$year)], na.rm = TRUE)
    rm(WHOdf)
    #Read in the Incidence spreadsheet and update
    IncidenceTable <- read.csv(file = paste0(datapath, '/incidence/VIMC-47-countries-for-cholera-modelling.csv'))
    NRCountryIndex <- match(country, IncidenceTable$country_code)
    
    YearList <- paste0(YearList[[1]], sep = '-')
    YearList[nyears] <- strsplit(YearList[nyears], '-')[[1]][1]
    for (i in 2:nyears){
      YearList[i] <- paste0(YearList[i-1], YearList[i])
    }
    
    IncidenceTable$year_list[NRCountryIndex] <- YearList[nyears]
    IncidenceTable$singular_estimate[NRCountryIndex] <- AverageAnnualCases
    IncidenceTable$is_num_case[NRCountryIndex] <- 1
    write.csv(IncidenceTable, file = paste0(datapath, '/incidence/VIMC-47-countries-for-cholera-modelling.csv'), row.names=FALSE)
  }
  
  #return IncidenceTable
  return(IncidenceTable)
}
