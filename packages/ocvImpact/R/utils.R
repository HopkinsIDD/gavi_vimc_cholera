
#' @name create_relative_worldpop_weights
#' @title create_relative_worldpop_weights
#' @description Generate a raster with the relative population weight of each cell, which will be used to extrapolate rasterized populations using the montagu population data
#' @param datapath path to data
#' @param country country code
#' @return raster of relative population weights according to worldpop
#' @export
create_relative_worldpop_weights <- function(datapath, country){

  pop <- load_worldpop_by_country(datapath, country)
  total_pop <- sum(raster::getValues(pop), na.rm = TRUE)

  wts <- raster::calc(pop, fun=function(x) x/total_pop)

  if(round(sum(raster::values(wts), na.rm = TRUE),10) != 1){
    warning(paste("Relative population weights do not sum to 1 in", country))
  }

  return(wts)
}


#' @name create_model_pop_raster
#' @title create_model_pop_raster
#' @description Generate a population raster from relative worldpop weights and standardized country population from Montagu
#' @param datapath path to input data
#' @param modelpath path to montagu files
#' @param country country code
#' @param year population year
#' @param cache montagu cache
#' @return raster of model population
#' @export
#' @include utils_montagu.R
create_model_pop_raster <- function(datapath, modelpath, country, year, cache){

  wts <- create_relative_worldpop_weights(datapath, country)
  year_pop <- import_country_population_1yr(modelpath, country, year, cache)

  rc <- raster::calc(wts, fun=function(x) x*year_pop)

  return(rc)
}


#' @name allocate_vaccine
#' @title allocate_vaccine
#' @description Create dataframe of total population vaccination coverage according to the admin-level assignments from assign_vaccine_targets
#' @param datapath path to input data
#' @param modelpath path to montagu files
#' @param country country code
#' @param scenario Unique string that identifies the coverage scenario name
#' @param cache montagu cache
#' @param ... Optional parameters to pass to [`assign_vaccine_targets()`]. See [`assign_vaccine_targets()`] for defaults.
#' @importFrom magrittr %>%
#' @return dataframe with proportion of total population allocated with vaccines in admin units in a given country and year
#' @export
#' @include utils_targeting.R load_shapefile_by_country.R utils.R
allocate_vaccine <- function(datapath, modelpath, country, scenario, cache, ...){

  vacc_targets <- assign_vaccine_targets(datapath, modelpath, country, scenario, cache, ...)

  ## skip if no vaccination
  if (is.null(vacc_targets)){
    vacc_coverage <- NULL
  } else{
    vacc_years <- sort(unique(vacc_targets$vacc_year))
    shp <- load_shapefile_by_country(datapath, country)

    ### a little play on the dataframe -- 7/2021
    shp <- shp %>%
      dplyr::mutate(genID = paste0(NAME_0, '-', NAME_1, '-', NAME_2))
    shp_sp <- GADMTools::gadm_sp_loadCountries(c(country), level = 2, basefile = file.path(datapath, "shapefiles/"))$spdf
    shp_sp$genID <- paste0(shp_sp$NAME_0, '-', shp_sp$NAME_1, '-', shp_sp$NAME_2)
    shp <- merge(shp, shp_sp, id = 'genID')

    vacc_pop <- lapply(vacc_years, function(yr){
      model_pop_raster <- create_model_pop_raster(datapath, modelpath, country, yr, cache)
      model_pop_admin <- get_admin_population(model_pop_raster, shp)
      rc <- shp %>%
        dplyr::mutate(vacc_year = yr,
                      pop_model = model_pop_admin) %>%
        dplyr::select(GID_2, vacc_year, pop_model)

      return(rc)
    }) %>%
      data.table::rbindlist()

    vacc_coverage <- dplyr::left_join(vacc_targets, vacc_pop, by = c("GID_2", "vacc_year")) %>%
      dplyr::mutate(actual_prop_ocv1_vaccinated = actual_ocv1_fvp/pop_model) %>%
      dplyr::mutate(actual_prop_ocv2_vaccinated = actual_ocv2_fvp/pop_model) %>%
      dplyr::mutate(actual_prop_atleast_1dose_vaccinated = (actual_ocv1_fvp + actual_ocv2_fvp)/pop_model) %>%
      sf::st_as_sf()

    if(any(vacc_coverage$actual_prop_ocv1_vaccinated>1 | vacc_coverage$actual_prop_ocv2_vaccinated>1 | vacc_coverage$actual_prop_atleast_1dose_vaccinated>1)){
      warning(paste("Number of fully vaccinated persons exceeds population in some admin units of", country))
    }
  }

  return(vacc_coverage)
}


#' @name get_model_years
#' @title get_model_years
#' @description Create list of years where vaccine dynamics will be important
#' @param modelpath path to montagu files
#' @param country country code
#' @param vacc_alloc object returned from [`allocate_vaccine()`]
#' @param cache montagu cache
#' @return list with model_years, start_year, and real_model_years
#' @export
#' @include utils_montagu.R
get_model_years <- function(modelpath, country, vacc_alloc, cache){
  tmp <- import_centralburden_template(mpathname, country, cache, redownload = FALSE)
  max_output_year <- max(tmp$year)

  if (!is.null(vacc_alloc)){
    ## assume that vaccine effects could be observed for maximum 8 years after the last campaign
    myear <- seq(min(as.numeric(vacc_alloc$vacc_year)), min(max(as.numeric(vacc_alloc$vacc_year)+10), max_output_year)) ## COULD CHANGE 10YEARS THIS TO LINK WITH MY_TRUNC_YEAR PARAM IN GENERATE_PCT_PROTECT_FUNCTION
  } else{
    myear <- NULL
  }

  myears_ls <- list(model_years = myear, output_years = sort(unique(tmp$year)))

  return(myears_ls)
}



#' @name generate_pct_protect_function
#' @title generate_pct_protect_function
#' @description Using vaccine efficacy studies in Bi et al. (2017), generate a function that takes the year since vaccination and provides an estimate of the direct ve. Vaccine efficacy declines to 0 after my_trunc_year years
#' @param my_trunc_year number of years after which VE is 0
#' @param my_ve_scen VE scenario
#' @importFrom magrittr %>%
#' @return dataframe with estimates for unweighted and weighted mean vaccination campaign coverage proportion across coverage surveys in the review
#' @export
generate_pct_protect_function <- function(my_trunc_year = 5, my_ve_scen = "base"){

    ve.dat <- readr::read_csv("input_data/ocv_ve_overtime.csv")

    ve.dat$T <- (ve.dat$TL+ve.dat$TR)/2

    ve.dat$weights <-  1/(abs(ve.dat$se)^2)

    ve.dat <- ve.dat[ve.dat$TL<48,]

    ##this is our basic trend in vaccine effictiveness.
    ve.trend <- lm(yi~T , data=ve.dat, weights = weights)

    ##create a function that gives the expected percent protected
    ##by a vaccine during a particular year after vaccination

    pct.protect<-function(year, ci=FALSE,trunc_year=my_trunc_year,ve_scen=my_ve_scen) {
        if (any(year%%1 !=0)) {warning("function designed to average across years only")}
        if(!ve_scen %in% c("base","low","high")) {warning("ve_scen not recognized, assuming ve_scen='base'")}

        months <- (year-.5)*12

        if(ve_scen %in% c("high","low")){
            ci <- TRUE
        }

        if (ci) {
            rc <- pmax(1-exp(predict(ve.trend, newdata=data.frame(T=months), interval="confidence") ),0) %>% unname
        } else {
            rc <- pmax(1-exp(predict(ve.trend, newdata=data.frame(T=months))),0) %>% unname
        }

        rc <- as.matrix(rc)

        if(any(year>trunc_year)){
            rc[which(year>trunc_year),] <- 0
        }

        if(ve_scen == "high"){
            rc <- rc[,2]
        } else if (ve_scen == "low"){
            rc <- rc[,3]
        }

        return(rc)
    }
    return(pct.protect)
}

#' @title pct_protect_all_ages_1d
#' @name pct_protect_all_ages_1d
#' @description 1 dose VE function for all age groups
#' @param years number of years since vacc campaign
#' @param proportion_under5 proportion of population under 5 years
#' @param my_trunc_year number of years after which VE is 0, default = 10
#' @return proportion protected for a given time since vaccination and proportion under 5
#' @export
pct_protect_all_ages_1d <- function(years, proportion_under5, my_trunc_year = 10){

  ve_under5_1d <- pct_protect_under5_1d(years, my_trunc_year)
  ve_over5_1d <- pct_protect_over5_1d(years)
  vaccine_efficacy <- (ve_under5_1d * proportion_under5 + ve_over5_1d * (1-proportion_under5))

  return(vaccine_efficacy)
}

#' @title pct_protect_under5_1d
#' @name pct_protect_under5_1d
#' @description 1 dose VE function for under 5s
#' @param years number of years since vacc campaign
#' @param my_trunc_year number of years after which VE is 0, default = 10
#' @return proportion protected for a given time since vaccination and under 5s
#' @export
pct_protect_under5_1d <- function(years, my_trunc_year = 10){
  ## currently assuming no protection for under 5 1-dose recipients
  return(0)
}

#' @title pct_protect_over5_1d
#' @name pct_protect_over5_1d
#' @description 1 dose VE function for over 5s
#' @param years number of years since vacc campaign
#' @param my_trunc_year number of years after which VE is 0, default = 10
#' @return proportion protected for a given time since vaccination and over 5s
#' @export
pct_protect_over5_1d <- function(years, my_trunc_year = 10){
  if (any(years%%1 !=0)) {warning("function designed to average across years only")}

  ve_1d_trend <- readRDS("input_data/log1d_obs.rds")
  months <- (years-.5)*12
  df <- metafor::predict.rma(ve_1d_trend, newmods=months, transf=function(x) 1-exp(x), addx=TRUE) %>%
    as.data.frame() %>%
    dplyr::mutate(pred_corr = ifelse(pred < 0, 0, pred))
  rc <- df$pred_corr

  if(any(years>my_trunc_year)){
    rc[which(years>my_trunc_year)] <- 0
  }

  return(rc)
}

#' @title pct_protect_all_ages_2d
#' @name pct_protect_all_ages_2d
#' @description 2 dose VE function for all age groups
#' @param years number of years since vacc campaign
#' @param proportion_under5 proportion of population under 5 years
#' @param my_trunc_year number of years after which VE is 0, default = 10
#' @return proportion protected for a given time since vaccination and proportion under 5
#' @export
pct_protect_all_ages_2d <- function(years, proportion_under5, my_trunc_year = 10){

  ve_under5_2d <- pct_protect_under5_2d(years, my_trunc_year)
  ve_over5_2d <- pct_protect_over5_2d(years, my_trunc_year)
  vaccine_efficacy <- (ve_under5_2d * proportion_under5 + ve_over5_2d * (1-proportion_under5))

  return(vaccine_efficacy)
}

#' @title pct_protect_under5_2d
#' @name pct_protect_under5_2d
#' @description 2 dose VE function for under 5s
#' @param years number of years since vacc campaign
#' @param my_trunc_year number of years after which VE is 0, default = 10
#' @return proportion protected for a given time since vaccination and under 5s
#' @export
pct_protect_under5_2d <- function(years, my_trunc_year = 10){
  age_multiplier <- mean(c(0.28, 1.04, 0.49, 0.83, 0.69, 0.56, 1.03, 0.34)) ## raw mean of relative VE of <5s to >5s in prelim version of Xu et al.
  rc <- age_multiplier * pct_protect_over5_2d(years, my_trunc_year)
  return(rc)
}

#' @title pct_protect_over5_2d
#' @name pct_protect_over5_2d
#' @description 2 dose VE function for over 5s
#' @param years number of years since vacc campaign
#' @param my_trunc_year number of years after which VE is 0, default = 10
#' @return proportion protected for a given time since vaccination and over 5s
#' @export
pct_protect_over5_2d <- function(years, my_trunc_year = 10){
  if (any(years%%1 !=0)) {warning("function designed to average across years only")}

  ve_2d_trend <- readRDS("input_data/log2d_obs.rds")
  months <- (years-.5)*12
  df <- metafor::predict.rma(ve_2d_trend, newmods=months, transf=function(x) 1-exp(x), addx=TRUE) %>%
    as.data.frame() %>%
    dplyr::mutate(pred_corr = ifelse(pred < 0, 0, pred))
  rc <- df$pred_corr

  if(any(years>my_trunc_year)){
    rc[which(years>my_trunc_year)] <- 0
  }

  return(rc)
}


#' @title pct_protect_all
#' @name pct_protect_all
#' @description total VE function for all age groups
#' @param years number of years since vacc campaign
#' @param proportion_under5 proportion of population under 5 years
#' @param proportion_one_dose proportion of one dose vaccinees
#' @param my_trunc_year number of years after which VE is 0, default = 10
#' @return proportion protected for a given time since vaccination and proportion under 5
#' @export
pct_protect_all <- function(years, proportion_under5, proportion_one_dose, my_trunc_year = 10){

  ve_2d <- pct_protect_all_ages_2d(years, proportion_under5, my_trunc_year)
  ve_1d <- pct_protect_all_ages_1d(years, proportion_under5, my_trunc_year)

  if(is.na(proportion_one_dose) | proportion_one_dose > 1 | proportion_one_dose < 0){
    message("proportion_one_dose is invalid, there may be no vaccination. returning VE = 0")
    vaccine_efficacy <- 0
  } else{
    vaccine_efficacy <- (ve_1d * proportion_one_dose + ve_2d * (1-proportion_one_dose))
  }

  return(vaccine_efficacy)
}

#' @name get_pop_proportion_ocv1
#' @title get_pop_proportion_ocv1
#' @description function that gets vacc_alloc and year as inputs and returns the proportion of vaccinated people that got
#' one dose of the vaccine
#' @param vacc_alloc the output of allocate_vaccine
#' @param year the vaccination year
#' @return the proportion of vaccinated people that received one dose of OCV in that year or the last year with vaccination
#' @export
get_pop_proportion_ocv1 <- function(vacc_alloc, year){
  vacc_alloc$vacc_year <- as.numeric(vacc_alloc$vacc_year)
  vacc_alloc_year <- vacc_alloc[vacc_alloc$vacc_year == year,]

  if(year < min(vacc_alloc$vacc_year)){
    warning(paste("Cannot return ocv1 proportion for years before the first vaccination year. Returning proportion for first vaccination year.", min(vacc_alloc$vacc_year)))
    vacc_alloc_year <- vacc_alloc[vacc_alloc$vacc_year == min(vacc_alloc$vacc_year),]

  } else if(nrow(vacc_alloc_year) == 0){
    # otherwise, if no vaccination in year, identify most recent vacc year prior to current year
    most_recent_vacc_year <- max(vacc_alloc$vacc_year[which(vacc_alloc$vacc_year < year)])
    vacc_alloc_year <- vacc_alloc[vacc_alloc$vacc_year == most_recent_vacc_year,]
    message(paste("There was no vaccination in", year, ". Calculating prop in", most_recent_vacc_year))
  }

  total_ocv1 <- sum(vacc_alloc_year$actual_ocv1_fvp) ## fvps with one dose
  total_ocv2 <- sum(vacc_alloc_year$actual_ocv2_fvp) ## fvps with two doses
  prop_ocv1_vs_ocv2 <- total_ocv1/(total_ocv1 + total_ocv2)

  return(prop_ocv1_vs_ocv2)
}



#' @name generate_indirect_incidence_mult
#' @title generate_indirect_incidence_mult
#' @description Using studies on indirect vaccine protection from Kolkota and Matlab, generate a function that takes the level of vaccination coverage and returns a multiplier indicating the percentage reduction in incidence due to indirect vaccine protection. These represent the baseline parameters for indirect vaccine protection.
#' @importFrom magrittr %>%
#' @return function for indirect protection
#' @export
generate_indirect_incidence_mult <- function(){

    ## data from Kolkata trial
    ## Herd protection by a bivalent killed whole-cell oral cholera vaccine in the slums of Kolkata, India
    vc_u_k <- c(25,28,31,34,45) # upper limit of coverage bins
    vc_l_k <- c(0,25,28,31,34) # lower limit of coverage bins
    risk_v_k <- c(1.28,1.48,1.01,0.97,1.24) # incid among vacc recipients
    risk_p_k <- c(5.54,5.64,2.48,2.25,1.93) # incid among placebo recipients

    ## data from Matlab trial
    ## also check this paper? Herd protection of unvaccinated adults by oral cholera vaccines in rural Bangladesh
    vc_u_m <- c(28,35,40,50,60) # upper limit of coverage bins
    vc_l_m <- c(0,28,35,40,50) # lower limit of coverage bins
    risk_v_m <- c(2.66,2.47,1.57,2.25,1.27) # incid among vacc recipients
    risk_p_m <- c(7.01,5.87,4.72,4.65,1.47) # incid among placebo recipients

                                        #effective_coverage <- function(cov){ cov*.75 + .25}

    df = data.frame(coverage=c((vc_u_k+vc_l_k)/2, (vc_u_m+vc_l_m)/2)/100,
                    vaccinated=c(risk_v_k, risk_v_m),
                    placebo=c(risk_p_k, risk_p_m),
                    ve_I=c(1-risk_v_k/risk_p_k, 1-risk_v_m/risk_p_m),
                    loc=c(rep("Kolkata", length(vc_u_k)), rep("Matlab", length(vc_u_m))))

    df <- dplyr::group_by(df, loc) %>%
      dplyr::mutate(indirect=pmax(0.001, 1 - (placebo/placebo[1])), effective_cov=coverage) %>%
      dplyr::ungroup()

    ##df %>%  ggplot(aes(x=coverage,y=indirect)) + geom_point(aes(color=loc))

    fit <- lm(log(indirect/(1-indirect)) ~ effective_cov, data=df)
    pred <- predict(fit, newdata = data.frame(effective_cov=seq(0.01,1,length=100)))

    indirect_inc_mult <- function(effective_coverage){
        tmp <- predict(fit, newdata = data.frame(effective_cov=effective_coverage))
        rc <- (1- (exp(tmp)/(1+exp(tmp)))) %>% unname

        rc[which(effective_coverage==0)] <- 1
        rc[which(effective_coverage==1)] <- 0
        # returns a multiplier for incidence based on coverage
        # (if effective coverage is 0%, there is no reduction due to indirect effects and multiplier is 1)
        return(rc)
    }

}


#' @name generate_flatline_multiplier
#' @title generate_flatline_multiplier
#' @description Generate a function that represents projected secular trends in cholera incidence. This multiplier is used to adjust projected cholera incidence in future years.
#' @param trendtype what type of trend function to return
#' @param datapath path to input data
#' @param modelpath path to montagu files
#' @param country country code
#' @param use_country_incid_trend whether use the country-level incidence rate trend
#' @return a function that can take in year or other parameters to produce a value
#' @export
#' @include incidence_rate_trend_multiplier.R
generate_flatline_multiplier <- function(trendtype,
                                         datapath,
                                         modelpath,
                                         country,
                                         use_country_incid_trend){
  #flatline_multiplier <- function(year, base_year = 2016) { return(1) }
  #return(flatline_multiplier)
  print('Now beginning using the generate_flatline_multiplier function. ')
  ### 11/03/2021 change -- now it works as a checkpoint
  if(trendtype == 'incidence rate'){
    return(ocvImpact::incidence_rate_trend_multiplier(datapath,
                                                      modelpath,
                                                      country,
                                                      use_country_incid_trend = use_country_incid_trend))

    # FunctionName <- paste0('incidence_rate_trend_multiplier_function_', country)
    # if(exists(FunctionName)){
    #   return(get(FunctionName))
    # }else{
    #   assign(paste0('incidence_rate_trend_multiplier_function_', country),
    #           ocvImpact::incidence_rate_trend_multiplier(datapath,
    #                                                      modelpath,
    #                                                      country,
    #                                                      use_country_incid_trend = use_country_incid_trend))
    #   return(get(FunctionName))

    # }
  }

}


#' @name generate_cfr
#' @title generate_cfr
#' @description Generate rough estimate of cholera case-fatality ratio based on previous WHO data
#' @param country country code
#' @importFrom magrittr %>%
#' @return numeric value of cfr (deaths/cases)
#' @export
generate_cfr <- function(country){

  deaths_summary <- readr::read_csv("input_data/who_cfrs.csv") %>%
    dplyr::select(-country_name)

  total_cfrs <- deaths_summary %>%
    dplyr::filter(!is.na(cases), !is.na(deaths), cases>0) %>%
    dplyr::filter(cfr <= 0.07)

  if (country %in% total_cfrs$cntry_code){
    calcs <- dplyr::filter(total_cfrs, cntry_code==country) %>%
      dplyr::summarise(cases = sum(cases), deaths = sum(deaths)) %>%
      dplyr::mutate(cfr = deaths/cases)

  } else {
    calcs <- dplyr::summarise(total_cfrs, cases = sum(cases), deaths = sum(deaths)) %>%
      dplyr::mutate(cfr = deaths/cases)
  }

  cfr <- calcs$cfr

  return(cfr)

}


#' @name generate_aoi
#' @title generate_aoi
#' @description Generate average age of infection
#' @param country country code
#' @return numeric value of cfr (deaths/cases)
#' @export
generate_aoi <- function(country){

    ## values for average age derived from this code:
    ## library(metafor)
    ## age_cases <- read_xlsx("data/age_cases.xlsx")

    ## ## using relationship between variacne and mean of age to get variance of log(age)
    ## ## 10.1002/sim.1525
    ## age_cases %>% filter(!is.na(age_sd)) %>%
    ## mutate(mean_log=log(mean_age),var_log=log(1+age_sd^2/mean_age^2)) %>%
    ## summarize(mean_age_inf=weighted.mean(mean_age,w=1/var_log)) %>% unlist

    ## aoi_south_asia <-  age_cases %>% filter(!is.na(age_sd),location %in% c('IND','BGD')) %>%
    ## mutate(mean_log=log(mean_age),var_log=log(1+age_sd^2/mean_age^2)) %>%
    ## summarize(mean_age_inf=weighted.mean(mean_age,w=1/var_log)) %>% unlist

    ## aoi_africa_haiti <- age_cases %>% filter(!is.na(age_sd),!location %in% c('IND','BGD')) %>%
    ## mutate(mean_log=log(mean_age),var_log=log(1+age_sd^2/mean_age^2)) %>%
    ## summarize(mean_age_inf=weighted.mean(mean_age,w=1/var_log)) %>% unlist


    aoi <- dplyr::case_when(
      country %in% c("BGD","IND","NPL") ~ 22.42,
      !country %in% c("BGD","IND","NPL") ~ 25.75
    )

    return(aoi)
}


#' @name generate_infectionDuration
#' @title generate_infectionDuration
#' @description Generate infection duration in years
#' @return numeric value of infection duration
#' @export
generate_infectionDuration <- function(){
  return(4/365)
}


#' @name load_custom_shapefile_by_country
#' @title load_custom_shapefile_by_country
#' @description Load a custom shapefile for DRC health zones (for the DRC case study) or download or admin0 (simple) level shapefile for a country from GADM, if not already in shapefiles folder, and return sf object
#' @param datapath path to data 
#' @param country country code
#' @param simple logical indicating whether simple (full country) shapefile should be used (default: FALSE)
#' @return shapefile in sf format
#' @export 
load_custom_shapefile_by_country <- function(datapath, country, simple = FALSE){
  
  
  tryCatch(
    {
      if (simple){
        country_pattern <- paste(country, "0", sep = "_")
        shp <- GADMTools::gadm_sf_loadCountries(c(country), level = 0, basefile = file.path(datapath, "shapefiles/"))$sf
        message(paste0("Loading ", datapath, "/shapefiles/", country_pattern, "_sf.rds"))
      } else{ 
        shp <- readRDS(config$shapefile_filename)
        message(paste0("Loading ", datapath, "/shapefiles/", config$shapefile_filename))
      }
    },
    error = function(e) {
      print(paste("Unable to get shapefile for", country, ".", e))
    }
  )
  
  return(shp)
}