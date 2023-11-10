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
#' @param clean logical that indicates whether existing targeting and model outputs (sus, pop, vacc) should be deleted (default = TRUE)
#' @param redraw logical that indicates whether existing incidence raster samples should be redrawn (default = FALSE)
#' @param ... Optional parameters to pass to [`assign_vaccine_targets()`]. See [`assign_vaccine_targets()`] for defaults.
#' @return 
#' @export
#' @include utils.R create_static_modelInputs.R create_sus_modelInputs.R create_expectedCases.R
run_country_scenario <- function(
  datapath, 
  modelpath, 
  country, 
  scenario,
  rawoutpath,
  nsamples, 
  num_doses = NULL,
  ve_direct = generate_pct_protect_function(),
  indirect_mult = generate_indirect_incidence_mult(),
  secular_trend_mult = function(a,b,c,d){return(a*b*c*d)},
  clean = TRUE,
  redraw = FALSE,
  ...){

  ##calam added conditionals changing the default ve_direct function to the new functions for the one and two-dose scenarios for the 2023 touchstone
  if (num_doses == "one" & scenario == "campaign-default"){ #use the new one dose vaccine efficacy function for the 2023 touchstones
    ve_direct = generate_pct_protect_function_one_dose()
  } else if (num_doses =="two" & scenario == "campaign-default"){ #use the new two dose vaccine efficacy function for the 2023 touchstones
    ve_direct = generate_pct_protect_function_two_dose()
  }
  ##finished adding conditionals, rest of the function is the same as the pre-2023 version
  
  vacc_alloc <- allocate_vaccine(datapath, modelpath, country, scenario, ...) #the changes start from here 

  ## write proportion vaccinated to file and export total population raster stack
  dummy <- create_static_modelInputs(datapath, modelpath, country, scenario, rawoutpath, vacc_alloc, clean)

  ## write susceptible population proportion raster 
  dummy2 <- create_sus_modelInputs(datapath, modelpath, country, scenario, rawoutpath, vacc_alloc, ve_direct, clean)

  incidence_rate_trend <- as.logical(config$setting$incidence_rate_trend)
  outbreak_multiplier <- as.logical(config$setting$outbreak_multiplier)
  setting <- paste0('incid_trend_', incidence_rate_trend, '_outb_layer_',  outbreak_multiplier)

  dir.create(paste0(rawoutpath, "/", scenario, "/", setting), showWarnings = FALSE)
  ec_out_fn <- paste0(rawoutpath, "/", scenario, "/", setting, "/", country, "_ec.csv")
  if(clean | !file.exists(ec_out_fn)){ ## rerun

    if (is.null(vacc_alloc)){
      message("Calculate expected cases: with vaccination")
      expCases <- create_expectedCases(datapath, modelpath, country, scenario, rawoutpath, vacc_alloc, indirect_mult, secular_trend_mult, nsamples, is_cf = TRUE, redraw)
    } else{
      message("Calculate expected cases: no vaccination")
      expCases <- create_expectedCases(datapath, modelpath, country, scenario, rawoutpath, vacc_alloc, indirect_mult, secular_trend_mult, nsamples, is_cf = FALSE, redraw)
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



