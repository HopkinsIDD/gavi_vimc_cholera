#' @name run_country_scenario
#' @title run_country_scenario
#' @description Run full country scenario
#' @param datapath path to input data 
#' @param modelpath path to montagu files
#' @param country country code
#' @param scenario Unique string that identifies the coverage scenario name
#' @param mu inverse life expectancy for country in years^(-1)
#' @param ve_direct vaccine effect function (default: generate_pct_protect_function())
#' @param indirect_mult indirect effects function (default: generate_indirect_incidence_mult())
#' @param secular_trend_mult function that returns a multiplier for secular incidence trends over time (default: generate_flatline_multiplier())
#' @param rawoutpath path to raw model output files
#' @param nsamples number of stochastic samples to use
#' @param clean logical that indicates whether existing vacc_files should be deleted (default = TRUE)
#' @param ... Optional parameters to pass to [`assign_vaccine_targets()`]. See [`assign_vaccine_targets()`] for defaults.
#' @return 
#' @export
run_country_scenario <- function(
  datapath, 
  modelpath, 
  country, 
  scenario,
  rawoutpath,
  nsamples, 
  ve_direct = generate_pct_protect_function(),
  indirect_mult = generate_indirect_incidence_mult(),
  secular_trend_mult = generate_flatline_multiplier(),
  clean = TRUE,
  ...){

  vacc_alloc <- allocate_vaccine(datapath, modelpath, country, scenario, ...)

  ## write proportion vaccinated to file and export total population raster stack
  dummy <- create_static_modelInputs(datapath, modelpath, country, scenario, rawoutpath, vacc_alloc, clean)

  ## write susceptible population proportion raster 
  dummy2 <- create_sus_modelInputs(datapath, modelpath, country, scenario, rawoutpath, vacc_alloc, ve_direct, clean)

  ec_out_fn <- paste0(rawoutpath, "/", scenario, "/", country, "_ec.csv")  
  if(clean | !file.exists(ec_out_fn)){ ## rerun

    if (is.null(vacc_alloc)){
      message("Calculate expected cases: with vaccination")
      expCases <- create_expectedCases(datapath, modelpath, country, scenario, rawoutpath, vacc_alloc, indirect_mult, secular_trend_mult, nsamples, is_cf = TRUE, clean)
    } else{
      message("Calculate expected cases: no vaccination")
      expCases <- create_expectedCases(datapath, modelpath, country, scenario, rawoutpath, vacc_alloc, indirect_mult, secular_trend_mult, nsamples, is_cf = FALSE, clean)
    } 

      ## Write to file 
      message(paste("Write expected cases:", country, scenario, "\n", ec_out_fn))
      readr::write_csv(expCases, ec_out_fn)

  } else{ ## read existing
    message(paste("Reading expected cases:", country, scenario, "\n", ec_out_fn))
    expCases <- readr::read_csv(ec_out_fn)
  }
  



  return(expCases)

}



