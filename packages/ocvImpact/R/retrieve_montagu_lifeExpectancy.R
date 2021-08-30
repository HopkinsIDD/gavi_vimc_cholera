#' @name retrieve_montagu_lifeExpectancy
#' @title retrieve_montagu_lifeExpectancy
#' @description Get life expectancy data from standardized source in Montagu.
#' @param modelpath path to Mongatu files
#' @param group_id the id of our modeling group (default = 'JHU-Lee')
#' @param expectations_all the expectations that VIMC wants for us (default = TRUE)
#' @param countries_all the countries that we want the total population data from (default = TRUE)
#' @return single or multiple CSV data files
#' @export 
retrieve_montagu_lifeExpectancy = function(modelpath, group_id = 'JHU-Lee', expectations_all = TRUE, countries_all = TRUE){
  
  
  ### Get the Touchstone from the Modelpath
  SplittedString = strsplit(modelpath, '/')[[1]]
  touchstone = SplittedString[length(SplittedString)]

  
  ### Commented out in 7/2021 since we have the montagu handle now
  ### Setup Montagu API (won't prompt anything if already logged in)
  #drat:::add("vimc")
  #montagu::montagu_server_global_default_set(
    #montagu::montagu_server("production", "montagu.vaccineimpact.org"))
  #invisible(montagu::montagu_scenarios(group_id, touchstone)) #this prompt is going to ask for username and password
  
  
  ##### Check the Expectations from Montagu and only keep the ones necessary -- modified in 7/2021
  ExpectationsIDList <- c()
  for (teams in montagu::montagu_expectations(group_id, touchstone)$description){
    if (group_id %in% strsplit(teams, ":")[[1]]){
      ExpectationsIDList <- c(ExpectationsIDList, montagu::montagu_expectations(group_id, touchstone)$id[match(teams, montagu::montagu_expectations(group_id, touchstone)$description)])
    }
  }
  
  
  ### Download the Life Expectancy Data
  if (expectations_all == TRUE){
    
    
    for (ExpectationsID in ExpectationsIDList) {
      ### Check the Countries Needed for Our Simulation
      CountriesList <- montagu::montagu_expectation_countries(group_id, touchstone, ExpectationsID)$id
      lx0 <- montagu::montagu_demographic_data("lx0", touchstone)
      
      if (countries_all == TRUE){
        FileName <- paste0(touchstone, '_', ExpectationsID, '_', 'all-countries', '_', 'lx0_both.csv')
        DirectoryFileName <- paste0(modelpath, '//', FileName)
        write.csv(lx0, DirectoryFileName, row.names = FALSE)
        rm(lx0) #to save the memory
      } else{
        lx0_subset <- subset(lx0, country_code %in% CountriesList)
        FileName <- paste0(touchstone, '_', ExpectationsID, '_', 'only-countries-needed', '_', 'lx0_both.csv')
        DirectoryFileName <- paste0(modelpath, '//', FileName)
        write.csv(lx0_subset, DirectoryFileName, row.names = FALSE)
        rm(lx0) #to save the memory
        rm(lx0_subset) #to save the memory
      }
      
    }

    
  } else{
    
    
    message(paste0("What life expectancy data with certain expectations that you want to download? Here are the available ones: "))
    print(montagu::montagu_expectations(group_id, touchstone))
    prompt <- "Please type the id's you want, with blank as seperation and no quotation marks: "
    SelectedExpectationsIDList <- unique(as.character(strsplit(readline(prompt), " ")[[1]]))
    
    #******if the input is wrong, ask again******
    while ( prod(as.integer(SelectedExpectationsIDList %in% as.character(ExpectationsIDList)))==0 & SelectedExpectationsIDList != 'cancel'){
      
      print(montagu::montagu_expectations(group_id, touchstone))
      prompt <- "Wrong input Format. Please type the id's you want, with blank as seperation between each one of them. No quotation mark is needed 
                (type <cancel> if you want to quit the whole simulation run, no progress will be saved): "
      SelectedExpectationsIDList <- unique(as.character(strsplit(readline(prompt), " ")[[1]]))
    }
    if (SelectedExpectationsIDList == 'cancel'){
      message(paste0("You quitted the whole run before you successfully retrieved the life expectancy data. The error message below is directly from the Montagu server: "))
    }
    
    SelectedExpectationsIDList = as.integer(SelectedExpectationsIDList)
    for (i in 1:length(SelectedExpectationsIDList)) {
      ExpectationsID <- SelectedExpectationsIDList[i]
      ### Check the Countries Needed for Our Simulation
      CountriesList <- montagu::montagu_expectation_countries(group_id, touchstone, ExpectationsID)$id
      lx0 <- montagu::montagu_demographic_data("lx0", touchstone)
      
      if (countries_all == TRUE){
        FileName <- paste0(touchstone, '_', ExpectationsID, '_', 'all-countries', '_', 'lx0_both.csv')
        DirectoryFileName <- paste0(modelpath, '//', FileName)
        write.csv(lx0, DirectoryFileName, row.names = FALSE)
        rm(lx0) #to save the memory
      } else{
        lx0_subset <- subset(lx0, country_code %in% CountriesList)
        FileName <- paste0(touchstone, '_', ExpectationsID, '_', 'only-countries-needed', '_', 'lx0_both.csv')
        DirectoryFileName <- paste0(modelpath, '//', FileName)
        write.csv(lx0_subset, DirectoryFileName, row.names = FALSE)
        rm(lx0) #to save the memory
        rm(lx0_subset) #to save the memory
      }
      
    }

    
  }
  message(paste0("The life expectancy data has been downloaded under ", modelpath))} 




