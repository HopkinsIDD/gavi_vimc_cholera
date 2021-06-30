#' @name retrieve_montagu_lifeExpectancy
#' @title retrieve_montagu_lifeExpectancy
#' @description Get life expectancy data from standardized source in Montagu.
#' @param modelpath path to Mongatu files
#' @param touchstone the unique id of the simulation
#' @param group_id the id of our modeling group (default = 'JHU-Lee')
#' @param expectations_all the expectations that VIMC wants for us (default = TRUE)
#' @param countries_all the countries that we want the total population data from (default = TRUE)
#' @return single or multiple CSV data files
#' @export 
retrieve_montagu_lifeExpectancy = function(modelpath, touchstone, group_id = 'JHU-Lee', expectations_all = TRUE, countries_all = TRUE){
  
  
  ### Check if the Files already Exist
  LifeExpFiles <- list.files(modelpath, pattern = "lx0_both.csv$")
  if (length(LifeExpFiles) >= 1){ #there should be 1 file
    
    return(message(paste0("The life expectancy data files have been under the directory: ", modelpath, '. No new download was made. ')))
  }
  
  
  ### Setup Montagu API
  drat:::add("vimc")
  montagu::montagu_server_global_default_set(
    montagu::montagu_server("production", "montagu.vaccineimpact.org"))
  invisible(montagu::montagu_scenarios(group_id, touchstone)) #this prompt is going to ask for username and password
  
  
  ##### Check the Expectations from Montagu
  ExpectationsIDList <- montagu::montagu_expectations(group_id, touchstone)$id
  
  
  ### Download the Life Expectancy Data
  if (expectations_all == TRUE){
    
    
    for (ExpectationsID in ExpectationsIDList) {
      ### Check the Countries Needed for Our Simulation
      CountriesList <- montagu::montagu_expectation_countries(group_id, touchstone, ExpectationsID)$id
      life_ex <- montagu::montagu_demographic_data("life_ex", touchstone)
      
      if (countries_all == TRUE){
        FileName <- paste0(touchstone, '_', ExpectationsID, '_', 'all-countries', '_', 'lx0_both.csv')
        DirectoryFileName <- paste0(modelpath, '\\', FileName)
        write.csv(life_ex, DirectoryFileName, row.names = TRUE)
        rm(life_ex) #to save the memory
      } else{
        life_ex_subset <- subset(life_ex, country_code %in% CountriesList)
        FileName <- paste0(touchstone, '_', ExpectationsID, '_', 'only-countries-needed', '_', 'lx0_both.csv')
        DirectoryFileName <- paste0(modelpath, '\\', FileName)
        write.csv(life_ex_subset, DirectoryFileName, row.names = TRUE)
        rm(life_ex) #to save the memory
        rm(life_ex_subset) #to save the memory
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
      life_ex <- montagu::montagu_demographic_data("life_ex", touchstone)
      
      if (countries_all == TRUE){
        FileName <- paste0(touchstone, '_', ExpectationsID, '_', 'all-countries', '_', 'lx0_both.csv')
        DirectoryFileName <- paste0(modelpath, '\\', FileName)
        write.csv(life_ex, DirectoryFileName, row.names = TRUE)
        rm(life_ex) #to save the memory
      } else{
        life_ex_subset <- subset(life_ex, country_code %in% CountriesList)
        FileName <- paste0(touchstone, '_', ExpectationsID, '_', 'only-countries-needed', '_', 'lx0_both.csv')
        DirectoryFileName <- paste0(modelpath, '\\', FileName)
        write.csv(life_ex_subset, DirectoryFileName, row.names = TRUE)
        rm(life_ex) #to save the memory
        rm(life_ex_subset) #to save the memory
      }
      
    }

    
  }
  message(paste0("The life expectancy data has been downloaded under ", modelpath))} 




