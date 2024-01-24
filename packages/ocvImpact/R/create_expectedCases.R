#' @name create_expectedCases
#' @title create_expectedCases
#' @description Create proportion of population vaccinated and total population rasterStacks. Write vaccination raster to file and export population raster
#' @param datapath path to input data 
#' @param modelpath path to montagu files
#' @param country country code
#' @param scenario Unique string that identifies the coverage scenario name
#' @param rawoutpath path to raw model output files
#' @param vacc_alloc object returned from [`allocate_vaccine()`]
#' @param indirect_mult function for calculating indirect effects
#' @param secular_trend_mult function for calculating secular trends in mean annual incidence
#' @param nsamples number of stochastic samples to use
#' @param is_cf logical indicating whether to run counterfactual model
#' @param redraw logical indicate whether to redraw incidence raster samples; pass to [`create_incid_raster()`]
#' @return dataframe with 
#' @export
#' @include utils.R utils_montagu.R create_incid_raster.R load_shapefile_by_country.R 
create_expectedCases <- function(
  datapath,
  modelpath,
  country,
  scenario,
  rawoutpath,
  vacc_alloc,
  indirect_mult,
  secular_trend_mult,
  nsamples,
  is_cf,
  redraw
  ){

  ############# Use the configs (might not work) -- 11/18/2021 #############
  incidence_rate_trend <- as.logical(config$setting$incidence_rate_trend)
  use_country_incid_trend <- as.logical(config$incid$use_country_incid_trend)
  outbreak_multiplier <- as.logical(config$setting$outbreak_multiplier)
  
  ## create template and inputs
  years_ls <- get_model_years(modelpath, country, vacc_alloc)
  model_years <- years_ls[["model_years"]]
  if (is_cf){
    model_years <- NULL
  }
  tmp <- import_centralburden_template(modelpath, country, redownload = FALSE)
  output_years <- sort(unique(tmp$year))

  sus_out_fn <- paste0(rawoutpath, "/", scenario, "/", country, "_sus.tif")
  vacc_out_fn <- paste0(rawoutpath, "/", scenario, "/", country, "_vacc.tif")
  pop_out_fn <- paste0(rawoutpath, "/", scenario, "/", country, "_pop.tif")

  ## write to file and import cholera incidence estimates
  lambda <- create_incid_raster(modelpath, datapath, country, nsamples, redraw)
  sus_rasterStack <- raster::brick(sus_out_fn)
  vacc_rasterStack <- raster::brick(vacc_out_fn)
  pop_rasterStack <- raster::brick(pop_out_fn)
  shp0 <- load_shapefile_by_country(datapath, country, simple = TRUE)
  
  #fixing the BGD issue
  if(country == 'BGD'){
    lambda <- raster::resample(lambda, pop_rasterStack, method = "ngb")
  }
  
  ### new tries -- 11/03/2021
  if(incidence_rate_trend){
    incid_trend_function <- ocvImpact::generate_flatline_multiplier(
                                        trendtype = 'incidence rate', 
                                        datapath = datapath, 
                                        modelpath = modelpath, 
                                        country = country, 
                                        use_country_incid_trend = use_country_incid_trend)
  }else{ 
    incid_trend_function <- function(year){return(1)}
  }
  
  # if(outbreak_multiplier & country != "NGA"){
  #   outbreak_trend_function <- ocvImpact::outbreak_incidence_rate_multiplier(
  #                                         datapath = datapath,
  #                                         modelpath = modelpath,
  #                                         country = country,
  #                                         lambda = lambda, 
  #                                         population_raster = pop_rasterStack, 
  #                                         output_years = output_years, 
  #                                         use_country_incid_trend = use_country_incid_trend)
      
  # }else{
  #   outbreak_trend_function <- function(yr_index){return(1)}
  # }
  # print('The outbreak_trend_function has been generated. ')

  ## make sure that the directory where all outbreak data is saved exists 
  dir.create(paste0(datapath, '/outbreak'), showWarnings = FALSE)

  ## apply outbreak multiplier only to the campaign years for now
  if (!scenario == "no-vaccination"){  ##using the following lines for the no-vaccination scenario creates inf values and leads to error
    first_year <- min(ocvImpact::import_coverage_scenario(modelpath, country, scenario, filter0 = FALSE, redownload = FALSE)$year)
    last_year  <- max(ocvImpact::import_coverage_scenario(modelpath, country, scenario, filter0 = FALSE, redownload = FALSE)$year) + 5
    campaign_years <- first_year:last_year  
  }

  #set.seed(666) #hopefully that the ramdom seed set by the previous steps could pass onto here and be used, but it needs to be checked
  ec_ls <- lapply(output_years, function(oy){

    ## calculate the outbreak multiplier for that year
    if(outbreak_multiplier){
      # first check if the outbreak raster file already exists, generate one if not
      random_seed <- as.numeric(config$setting$random_seed)
      setting_num <- random_seed #for now, yes
      outbreak_out_fn <- paste0(datapath, '/outbreak/', country, '_', (oy %% 10), '.tif')

      if(!file.exists(outbreak_out_fn)){
        outbreak_trend_function <- ocvImpact::outbreak_incidence_rate_multiplier(
                                              datapath = datapath,
                                              modelpath = modelpath,
                                              country = country,
                                              lambda = lambda, 
                                              population_raster = pop_rasterStack, 
                                              output_years = c(oy), 
                                              use_country_incid_trend = use_country_incid_trend)
      }else{
        outbreak_multiplier_raster <- raster::stack(outbreak_out_fn) #please note that the raster here is RasterStack, not RasterBrick
        outbreak_trend_function <- function(yr_index){return(outbreak_multiplier_raster)}
      }
      
    }else{
      outbreak_trend_function <- function(yr_index){return(1)}
    }



    yr_index <- which(oy == output_years)
    pop_rasterLayer <- raster::subset(pop_rasterStack, yr_index, drop = FALSE)
    vacc_rasterLayer <- raster::subset(vacc_rasterStack, yr_index, drop = FALSE)
    sus_rasterLayer <- raster::subset(sus_rasterStack, yr_index, drop = FALSE)

    incid_trend_multiplier <- incid_trend_function(year = as.numeric(output_years[yr_index]))
    outbreak_ic_multiplier <- outbreak_trend_function(yr_index = 1) ### the returned list only has one element to use anyways
    confirmation_multiplier <- 1 ### 11/25 added
    utilization_multiplier <- 1 ### 11/25 added
    
    overall_multiplier <- secular_trend_mult(incid_trend_multiplier, outbreak_ic_multiplier, confirmation_multiplier, utilization_multiplier)
    

    
    if (oy %in% model_years){ ## consecutive years where vaccine dynamics are in play

      ## make new indirect effects template
      indirect_rasterLayer <- pop_rasterLayer
      raster::values(indirect_rasterLayer) <- indirect_mult(1-as.numeric(raster::values(sus_rasterLayer)))

      ec_rasterStack <- tryCatch(
        if(!is.numeric(overall_multiplier) & class(overall_multiplier) == 'raster'){
          raster::overlay(
            sus_rasterLayer,
            pop_rasterLayer,
            lambda,
            indirect_rasterLayer,
            overall_multiplier, 
            fun = function(x, y, z, a, b){
              x*y*z*a*b
            },
            recycle = TRUE, unstack = TRUE) 
            
        }else {
          lambda * sus_rasterLayer * pop_rasterLayer  * indirect_rasterLayer * overall_multiplier
        },

        error = function(e){
          print(paste0('The year when it goes wrong is ', oy))
          print('This is a year when vaccination campaign is going on. ')
          print('The following is the class for each object: sus_rasterLayer, pop_rasterLayer, lambda, indirect_rasterLayer, and overall_multiplier: ')
          print(class(sus_rasterLayer))
          print(class(pop_rasterLayer))
          print(class(lambda))
          print(class(indirect_rasterLayer))
          print(class(overall_multiplier))
          print('The following is the nlayers for each object: sus_rasterLayer, pop_rasterLayer, lambda, and indirect_rasterLayer: ')
          print(raster::nlayers(sus_rasterLayer))
          print(raster::nlayers(pop_rasterLayer))
          print(raster::nlayers(lambda))
          print(raster::nlayers(indirect_rasterLayer))
          
        }
      )

    } else {
      
      ec_rasterStack <- tryCatch(
        if(!is.numeric(overall_multiplier) & class(overall_multiplier) == 'raster'){
          raster::overlay(
            pop_rasterLayer,
            lambda,
            overall_multiplier, 
            fun = function(x, y, z){
              return(x*y*z)
            },
            recycle = TRUE, unstack = TRUE)

        }else{
          lambda * pop_rasterLayer * overall_multiplier
        },

        error = function(e){
          print(paste0('The year when it goes wrong is ', oy))
          print('This is a year when no vaccination campaign is going on. ')
          print('The following is the class for each object: pop_rasterLayer, lambda, and overall_multiplier: ')
          print(class(pop_rasterLayer))
          print(class(lambda))
          print(class(overall_multiplier))
          print('The following is the nlayers/value for each object: pop_rasterLayer, lambda, and overall_multiplier: ')
          if(!is.numeric(overall_multiplier) & class(overall_multiplier) == 'raster'){
            print(list(lambda = raster::nlayers(lambda),
                       pop = raster::nlayers(pop_rasterLayer), 
                       overall_multiplier = raster::nlayers(overall_multiplier)))
          }else{
            print(list(lambda = raster::extent(lambda),
                       pop = raster::extent(pop_rasterLayer), 
                       overall_multiplier = overall_multiplier))
          }
          
        }
      )

    }
    print('ec_rasterStack has been generated')



    ### optimize memory usage
    rm(outbreak_ic_multiplier)
    rm(outbreak_trend_function)
    rm(overall_multiplier)
    gc()


    
    ### fixing the MRT issue
    if(country == 'MRT'){
      lambda <- raster::setExtent(lambda, raster::extent(shp0), keepres=FALSE, snap=FALSE)
      pop_rasterLayer <- raster::setExtent(pop_rasterLayer, raster::extent(shp0), keepres=FALSE, snap=FALSE)
    }
    


    ec_yr <- exactextractr::exact_extract(ec_rasterStack, shp0, fun = "sum", stack_apply = TRUE)
    print('The first exact_extract function got passed. ')

    mean_incid <- exactextractr::exact_extract(
      lambda, 
      shp0, 
      function(values, coverage_frac, weights){
        weighted.mean(values, ifelse(is.na(coverage_frac*weights), 0, coverage_frac*weights), na.rm = TRUE)
      },
      weights = pop_rasterLayer, 
      stack_apply = TRUE)
    print('The second exact_extract function got passed. ')

    ec_vec <- as.numeric(ec_yr)
    mean_incid_vec <- as.numeric(mean_incid)
    rc <- tibble::tibble(country = country, year = oy, run_id = seq_along(ec_vec), ec = ec_vec, incid_rate = mean_incid_vec)

    return(rc)

  })

  ec_final <- data.table::rbindlist(ec_ls)
  
  rm(lambda, sus_rasterStack, vacc_rasterStack, pop_rasterStack, ec_ls)
  gc()

  return(ec_final)
}