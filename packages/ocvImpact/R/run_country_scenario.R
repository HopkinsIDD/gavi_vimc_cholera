#' @name run_country_scenario
#' @title run_country_scenario
#' @description Run full country scenario
#' @param datapath path to input data
#' @param modelpath path to montagu files
#' @param country country code
#' @param scenario Unique string that identifies the coverage scenario name
#' @param mu inverse life expectancy for country in years^(-1)
#' @param indirect_mult indirect effects function (default: generate_indirect_incidence_mult())
#' @param secular_trend_mult function that returns a multiplier for secular incidence trends over time (default: generate_flatline_multiplier())
#' @param rawoutpath path to raw model output files
#' @param nsamples number of stochastic samples to use
#' @param num_doses the number of doses for the one-dose and two-dose scenarios of the 202310gavi-4 touchstone
#' @param clean logical that indicates whether existing targeting and model outputs (sus, pop, vacc) should be deleted (default = TRUE)
#' @param redraw logical that indicates whether existing incidence raster samples should be redrawn (default = FALSE)
#' @param ... Optional parameters to pass to [`assign_vaccine_targets()`]. See [`assign_vaccine_targets()`] for defaults.
#' @return dataframe and export files
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
    indirect_mult = generate_indirect_incidence_mult(),
    secular_trend_mult = function(a,b,c,d){return(a*b*c*d)},
    clean = TRUE,
    redraw = FALSE,
    ...){

  incidence_rate_trend <- as.logical(config$setting$incidence_rate_trend)
  outbreak_multiplier <- as.logical(config$setting$outbreak_multiplier)
  setting <- paste0('incid_trend_', incidence_rate_trend, '_outb_layer_',  outbreak_multiplier)
  dir.create(paste0(rawoutpath, "/", scenario, "/", setting), showWarnings = FALSE)
  
  ##calam added to write expected cases in separate directories for the one dose and two dose campaigns for 202310gavi-4 touchstone
  runname <- config$runname
  if (runname == "202310gavi-4"){
    ec_out_fn <- paste0(rawoutpath, "/", scenario, "/", setting, "/", country,"_",num_doses,"_ec.csv")
    cov_out_fn <- paste0(rawoutpath, "/", scenario, "/", setting, "/", country,"_",num_doses,"_coverage.csv") ##to write modelled coverage
    imm_out_fn <- paste0(rawoutpath, "/", scenario, "/", setting, "/", country,"_",num_doses,"_immune.csv") ##ca debug
  } else {
    ec_out_fn <- paste0(rawoutpath, "/", scenario, "/", setting, "/", country, "_ec.csv")
    cov_out_fn <- paste0(rawoutpath, "/", scenario, "/", setting, "/", country, "_coverage.csv")
    imm_out_fn <- paste0(rawoutpath, "/", scenario, "/", setting, "/", country, "_immune.csv")
  }
  if(clean | !file.exists(ec_out_fn)){ ## rerun
  
  
  ## avoid reading montagu files multiple times
  montagu_cache <- new.env()
  montagu_cache[["coverage_scenario"]] <- import_coverage_scenario(modelpath, country, scenario, montagu_cache, filter0 = FALSE, redownload = FALSE)
  montagu_cache[["centralburden_template"]] <- import_centralburden_template(modelpath, country, montagu_cache, redownload = FALSE)
  montagu_cache[["country_population"]] <- import_country_population(modelpath, country, montagu_cache, redownload = FALSE)
  montagu_cache[["int_country_population"]] <- import_int_country_population(modelpath, country, montagu_cache, redownload = FALSE)
  montagu_cache[["country_agePop"]] <- import_country_agePop(modelpath, country, montagu_cache, redownload = FALSE)
  montagu_cache[["country_lifeExpectancy"]] <- import_country_lifeExpectancy(modelpath, country, montagu_cache, redownload = FALSE)

  if (config$vacc$targeting == "custom"){
    ## use the custom-made targeting table for the 'custom' targeting strategy
    vacc_alloc <- readr::read_csv(paste0(datapath,"/custom_targeting.csv")) ##if the custom targeting table is inside /input_data
  } else {
    vacc_alloc <- allocate_vaccine(datapath, modelpath, country, scenario, montagu_cache, ...) #the changes start from here
  }

  ## write proportion vaccinated to file and export total population raster stack
  dummy <- create_static_modelInputs(datapath, modelpath, country, scenario, rawoutpath, vacc_alloc, montagu_cache, clean)

  ## write susceptible population proportion raster
  track_prop_immune <- create_sus_modelInputs(datapath, modelpath, country, scenario, rawoutpath, vacc_alloc, montagu_cache, clean)


    if (is.null(vacc_alloc)){
      message("Calculate expected cases: no vaccination")
      expCases <- create_expectedCases(datapath, modelpath, country, scenario, rawoutpath, vacc_alloc, indirect_mult, secular_trend_mult, nsamples, montagu_cache, is_cf = TRUE, redraw)
    } else{
      message("Calculate expected cases: with vaccination")
      expCases <- create_expectedCases(datapath, modelpath, country, scenario, rawoutpath, vacc_alloc, indirect_mult, secular_trend_mult, nsamples, montagu_cache, is_cf = FALSE, redraw)
    }

    ## Write to file
    message(paste("Write expected cases:", country, scenario, "\n", ec_out_fn))
    readr::write_csv(expCases, ec_out_fn)

    ##only write coverage for campaign scenarios (vacc_alloc is null for the no-vaccination scenarios)
    if (scenario == "ocv1-default" | scenario == "ocv1-ocv2-default"){
      message(paste("Write modelled coverage and prop immune:", country, scenario, "\n", cov_out_fn))
      readr::write_csv(vacc_alloc, cov_out_fn)
      readr::write_csv(track_prop_immune, imm_out_fn)
    }


  } else{ ## read existing
    message(paste("Reading expected cases:", country, scenario, "\n", ec_out_fn))
    expCases <- readr::read_csv(ec_out_fn)
  }




  return(expCases)

}



