#' @name retrieve_montagu_centralburden_template
#' @title retrieve_montagu_centralburden_template
#' @description Retrieves the CSV file(s) with central burden estimates template from the Montagu API.
#' @param modelpath path to Mongatu files
#' @param group_id the id of our modeling group (default = 'JHU-Lee')
#' @param expectations_all the expectations that VIMC wants for us (default = TRUE)
#' @return single or multiple CSV data files
#' @export 
retrieve_montagu_centralburden_template = function(modelpath, group_id = 'JHU-Lee', expectations_all = TRUE){
  
  
  ### Get the Touchstone from the Modelpath
  SplittedString = strsplit(modelpath, '/')[[1]]
  touchstone = SplittedString[length(SplittedString)]
  
  
  ### Setup Montagu API (won't prompt anything if already logged in)
  drat:::add("vimc")
  montagu::montagu_server_global_default_set(
    montagu::montagu_server("production", "montagu.vaccineimpact.org"))
  invisible(montagu::montagu_scenarios(group_id, touchstone)) #this prompt is going to ask for username and password
  
  
  ### Check the Expectations from Montagu
  ExpectationsIDList <- montagu::montagu_expectations(group_id, touchstone)$id
  
  
  ### Download the Template
  if (expectations_all == TRUE){
    
    for (ExpectationsID in ExpectationsIDList) {
      CentralBurdenTemplate <- montagu::montagu_central_burden_estimate_template(group_id, touchstone, ExpectationsID)
      #'central-burden-template.201910gavi-5.Cholera_ standard template.csv'
      FileName <- paste0('central-burden-template.', touchstone, '.Cholera_ standard template.', ExpectationsID, '.csv')
      DirectoryFileName <- paste0(modelpath, '\\', FileName)
      write.csv(CentralBurdenTemplate, DirectoryFileName, row.names = TRUE)
    }
    rm(CentralBurdenTemplate) #to save memory
    message(paste0("The central burden template has been downloaded under ", modelpath))
    
  } else{
    
    message(paste0("What specific template with certain expectations that you want to download? Here are the available ones: "))
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
      message(paste0("You quitted the whole run before you successfully retrieved the central burden template(s). The error message below is directly from the Montagu server: "))
    }
    
    
    SelectedExpectationsIDList = as.integer(SelectedExpectationsIDList)
    for (i in 1:length(SelectedExpectationsIDList)) {
      ExpectationsID <- SelectedExpectationsIDList[i]
      CentralBurdenTemplate <- montagu::montagu_central_burden_estimate_template(group_id, touchstone, ExpectationsID)
      FileName <- paste0('central-burden-template.', touchstone, '.Cholera_ standard template.', ExpectationsID, '.csv')
      DirectoryFileName <- paste0(modelpath, '\\', FileName)
      write.csv(CentralBurdenTemplate, DirectoryFileName, row.names = TRUE)
    }
    rm(CentralBurdenTemplate) #to save memory
    message(paste0("The central burden template has been downloaded under ", modelpath))
  }
}







