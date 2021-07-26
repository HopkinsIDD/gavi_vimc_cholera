#' @name retrieve_montagu_coverage
#' @title retrieve_montagu_coverage
#' @description Retrieves the CSV file(s) with vaccination coverage values for the scenario from the Montagu API.
#' @param modelpath path to Mongatu files
#' @param group_id the id of our modeling group (default = 'JHU-Lee')
#' @param scenarios_all the scenario that we want the coverage data from (default = TRUE)
#' @return single or multiple CSV data files
#' @export 
retrieve_montagu_coverage = function(modelpath, group_id = 'JHU-Lee', scenarios_all = TRUE){
  
  
  ### Get the Touchstone from the Modelpath
  SplittedString = strsplit(modelpath, '/')[[1]]
  touchstone = SplittedString[length(SplittedString)]
  
  
  ### Setup Montagu API (won't prompt anything if already logged in)
  drat:::add("vimc")
  montagu::montagu_server_global_default_set(
    montagu::montagu_server("production", "montagu.vaccineimpact.org"))
  invisible(montagu::montagu_scenarios(group_id, touchstone)) #this prompt is going to ask for username and password
  
  
  ### Check the Scenarios Available
  ScenariosList <- montagu::montagu_scenarios(group_id, touchstone)$scenario_id
  
  
  ### Download the Coverage Data
  if (scenarios_all == TRUE){
    
    for (Scenarios in ScenariosList) {
      cov <- montagu::montagu_coverage_data(group_id, touchstone, Scenarios)
      FileName <- paste0('coverage_', touchstone, '_', Scenarios, '.csv')
      DirectoryFileName <- paste0(modelpath, '//', FileName)
      write.csv(cov, DirectoryFileName, row.names = TRUE)
    }
    message(paste0("The coverage data has been downloaded under ", modelpath))
    
  } else{
    
    message(paste0("What specific scenarios you want to download? Here are the available ones: "))
    print(montagu::montagu_scenarios(group_id, touchstone))
    prompt <- "Please type the scenario_id's you want, with blank as seperation and no quotation marks: "
    SelectedScenariosList <- unique(as.character(strsplit(readline(prompt), " ")[[1]]))
    
    
    #******if the input is wrong, ask again******
    while ( prod(as.integer(SelectedScenariosList %in% ScenariosList))==0 & SelectedScenariosList != 'cancel'){
      
      print(montagu::montagu_scenarios(group_id, touchstone))
      prompt <- "Wrong input Format. Please type the scenario_id's you want, with blank as seperation between each one of them. No quotation mark is needed 
                (type <cancel> if you want to quit the whole simulation run, no progress will be saved): "
      SelectedScenariosList <- unique(as.character(strsplit(readline(prompt), " ")[[1]]))
    }
    if (SelectedExpectationsIDList == 'cancel'){
      message(paste0("You quitted the whole run before you successfully retrieved the coverage data. The error message below is directly from the Montagu server: "))
    }
    
    
    for (i in 1:length(SelectedScenariosList)) {
      Scenarios <- SelectedScenariosList[i]
      cov <- montagu::montagu_coverage_data(group_id, touchstone, Scenarios)
      FileName <- paste0('coverage_', touchstone, '_', Scenarios, '.csv')
      DirectoryFileName <- paste0(modelpath, '//', FileName)
      write.csv(cov, DirectoryFileName, row.names = TRUE)
    }
    message(paste0("The coverage data has been downloaded under ", modelpath))
  }
}







