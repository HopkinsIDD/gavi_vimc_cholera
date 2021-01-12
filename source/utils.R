
#' @name create_relative_worldpop_weights
#' @title create_relative_worldpop_weights
#' @description Generate a raster with the relative population weight of each cell, which will be used to extrapolate rasterized populations using the montagu population data
#' @param datapath path to data 
#' @param country country code
#' @return raster of relative population weights according to worldpop 
#' @export
create_relative_worldpop_weights <- function(datapath, country){

  pop <- load_worldpop_by_country(datapath, country)
  shp <- load_shapefile_by_country(datapath, country)
  pop_crop <- raster::crop(pop, shp, snap = "out")
  total_pop <- sum(raster::getValues(pop_crop), na.rm = TRUE)

  wts <- raster::calc(pop_crop, fun=function(x) x/total_pop)
  
  if(sum(raster::values(wts), na.rm = TRUE) != 1){
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
#' @return raster of model population
#' @export
create_model_pop_raster <- function(datapath, modelpath, country, year){

  wts <- create_relative_worldpop_weights(datapath, country)
  year_pop <- import_country_population_1yr(modelpath, country, year)

  rc <- raster::calc(wts, fun=function(x) x*year_pop)

  return(rc)
}


#' @name allocate_vaccine
#' @title allocate_vaccine
#' @description Create dataframe of total population vaccination coverage according to the admin-level assignments from assign_vaccine_targets
#' @param datapath path to input data 
#' @param modelpath path to montagu files
#' @param country country code
#' @param ... Optional parameters to pass to [`assign_vaccine_targets()`]. See [`assign_vaccine_targets()`] for defaults.
#' @importFrom magrittr %>%
#' @return dataframe with proportion of total population allocated with vaccines in admin units in a given country and year
#' @export
allocate_vaccine <- function(datapath, modelpath, country, ...){

  vacc_targets <- assign_vaccine_targets(datapath, modelpath, country, ...)

  vacc_years <- sort(unique(vacc_targets$vacc_year))
  shp <- load_shapefile_by_country(datapath, country)

  vacc_pop <- lapply(vacc_years, function(yr){
    model_pop_raster <- create_model_pop_raster(datapath, modelpath, country, yr)
    model_pop_admin <- get_admin_population(model_pop_raster, shp)
    rc <- sf::st_drop_geometry(shp) %>%
      dplyr::mutate(vacc_year = yr,
                    pop_model = model_pop_admin) %>%
      dplyr::select(GID_2, vacc_year, pop_model)
    
    return(rc)
  }) %>%
    data.table::rbindlist()

  vacc_coverage <- dplyr::left_join(vacc_targets, vacc_pop, by = c("GID_2", "vacc_year")) %>%
    dplyr::mutate(actual_prop_vaccinated = actual_fvp/pop_model)

  if(any(vacc_coverage$actual_vacc_prop>1)){
    warning(paste("Number of fully vaccinated persons exceeds population in some admin units of", country))
  }

  return(vacc_coverage)
}


#' @name generate_indirect_incidence_mult
#' @title generate_indirect_incidence_mult
#' @description Using studies on indirect vaccine protection from Kolkota and Matlab, generate a function that takes the level of vaccination coverage and returns a multiplier indicating the percentage reduction in incidence due to indirect vaccine protection. These represent the baseline parameters for indirect vaccine protection. 
#' @return 
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
      ungroup

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
#' @return
generate_flatline_multiplier <- function(){
  flatline_multiplier <- function(year, base_year = 2016) { return(1) }
  return(flatline_multiplier)
}
