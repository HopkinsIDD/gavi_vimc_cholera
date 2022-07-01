#' @name get_inc_tp_doses
#' @title get_inc_tp_doses
#' @description function to get incidence, target population and number of doses administered under selected campaign scenario, surveillance scenario and threshold
#' @param rawoutpath path of output_raw in which all the model output files are stored
#' @param configpath path of configuration files
#' @param countries vector of all the countries included the diagnostic report
#' @param scenario takes value of either "campaign-default" or "no-vaccination"
#' @param threshold threshold chosen for the diagnostic report, takes one of 1/1000, 1/5000, or 1/10,000
#' @param surveillance_scenario scale of surveillance, takes "no-estimate", "global-estimate", or "district-estimate"
#' @return 
#' @export
#' @include
get_inc_tp_doses <- function(rawoutpath = "output_raw/202110gavi-3",
                             configpath = "configs/202110gavi-3",
                             countries,
                             scenario,
                             surveillance_scenario,
                             threshold,
                             ...){
    
    # first, check whether desired output files are all in the output_raw folder 
    outputs_availability <- check_outputs_availability(rawoutpath = rawoutpath ,configpath = configpath, countries = countries,
                               scenario = scenario, surveillance_scenario = surveillance_scenario, threshold = threshold)
    if(outputs_availability == FALSE){
        stop(paste("Stop calculating target population as not all the required output files are available in", rawoutpath, "folder."))
    }else{

        # set up file path
        output_path <- paste0(rawoutpath, "/", scenario)
        df_setting <- check_model_setting(configpath = configpath, countries = countries,
                                          scenario = scenario, surveillance_scenario = surveillance_scenario, 
                                          threshold = threshold)
        
        # loop through countries
        for(country in countries){

            # get filenames for both admin levels
            incid <- df_setting[which(df_setting$country == country), ]$incidence_trend 
            outbk <- df_setting[which(df_setting$country == country), ]$outbreak_multiplier 
            vax_cov <- df_setting[which(df_setting$country == country), ]$vac_coverage
            surveillance_scenario <- df_setting[which(df_setting$country == country), ]$surveillance_scenario
            prefix <- paste0("incid_", incid, "_outbk_", outbk, "_")
            filename_admin1 <- paste0(output_path, "/", prefix, threshold, "_", surveillance_scenario, "_", country, "_rc_admin1.csv")
            filename_admin2 <- paste0(output_path, "/", prefix, threshold, "_", surveillance_scenario, "_", country, "_rc_admin2.csv")

            # read in csv files
            df_admin1 <- read.csv(filename_admin1)
            df_admin2 <- read.csv(filename_admin2)

            # if this is the first country, make new tables of target pop and doses administered
            if(which(countries == country) == 1){
                # when doing admin1 level OCV targeting
                df_admin1_sumadmins <- get_inc_tp_doses_helper(rc = df_admin1, sum_level = "country", vax_cov = vax_cov, surveillance_scenario = surveillance_scenario) # sum all admins together
                df_admin1_byadmin <- get_inc_tp_doses_helper(rc = df_admin1, sum_level = "admin", vax_cov = vax_cov, surveillance_scenario = surveillance_scenario) # by admin
                # when doing admin2 level OCV targeting
                df_admin2_sumadmins <- get_inc_tp_doses_helper(rc = df_admin2, sum_level = "country", vax_cov = vax_cov, surveillance_scenario = surveillance_scenario) # sum all admins together
                df_admin2_byadmin <- get_inc_tp_doses_helper(rc = df_admin2, sum_level = "admin", vax_cov = vax_cov, surveillance_scenario = surveillance_scenario) # by admin
            }else{ # else, append to existing result tables
                df_admin1_sumadmins <- rbind(df_admin1_sumadmins, get_inc_tp_doses_helper(rc = df_admin1, sum_level = "country", vax_cov = vax_cov, surveillance_scenario = surveillance_scenario))
                df_admin1_byadmin <- rbind(df_admin1_byadmin, get_inc_tp_doses_helper(rc = df_admin1, sum_level = "admin", vax_cov = vax_cov, surveillance_scenario = surveillance_scenario))
                df_admin2_sumadmins <- rbind(df_admin2_sumadmins, get_inc_tp_doses_helper(rc = df_admin2, sum_level = "country", vax_cov = vax_cov, surveillance_scenario = surveillance_scenario))
                df_admin2_byadmin <- rbind(df_admin2_byadmin, get_inc_tp_doses_helper(rc = df_admin2, sum_level = "admin", vax_cov = vax_cov, surveillance_scenario = surveillance_scenario))
            }
        }
    }
    # make a list storing all 4 result tables
    list_inc_tp_doses <- list("df_admin1_sumadmins" = df_admin1_sumadmins,
                              "df_admin1_byadmin" = df_admin1_byadmin,
                              "df_admin2_sumadmins" = df_admin2_sumadmins,
                              "df_admin2_byadmin" = df_admin2_byadmin)

    return(list_inc_tp_doses)
}

# test get_inc_tp_doses()
# a <- get_inc_tp_doses(countries = c("NGA"), scenario = "campaign-default", surveillance_scenario = "district-estimate", threshold = 0.001)




#' @name get_inc_tp_doses_helper
#' @title get_inc_tp_doses_helper
#' @description Function to assist calculation of get_inc_tp_doses(), read in only one data frame a time
#' @param rc a data frame of one scenario of one country
#' @param sum_level takes either "country" or "admin". If takes "admin", returns a table by admin units; if takes "country", returns a table for the whole country aggregating all admin units.
#' @param vax_cov vaccine coverage, default value is 0.68
#' @param surveillance_scenario scale of surveillance scenario, which will be used to calculate observed incidence 
#' @return a table with incidence, target population and number of doses administered, by year and cumulatively
#' @export
#' @include
get_inc_tp_doses_helper <- function(rc,
                                    sum_level,
                                    vax_cov = 0.68,
                                    surveillance_scenario){
  
  country <- rc$ISO[1]
  vax_years <- unique((rc %>% filter(!is.na(confirmation_lens)))$year)
  
  # fill in missing observed incidence column
  if(surveillance_scenario == "district-estimate"){
    rc <- rc %>% mutate(observed_incidence_rate = true_incidence_rate)
  }
  if(surveillance_scenario == "global-estimate"){
    global_comfirmation_rate <- unique(rc$confirmation_rate)[!is.na(unique(rc$confirmation_rate))]
    rc <- rc %>% mutate(observed_incidence_rate = true_incidence_rate / global_comfirmation_rate)
  }
  if(surveillance_scenario == "no-estimate"){
    rc <- rc %>% mutate(observed_incidence_rate = true_incidence_rate / confirmation_rate)
  }
  

  ### if want to summarize all admins together 
  if(sum_level == "country"){
    
    # by year
    df_tp_byyear <- rc %>% 
      filter(is_target == 1) %>%
      group_by(year, run_id) %>%
      summarize(target_pop = sum(pop_model)) %>%
      group_by(year) %>%
      summarize(tp_lb = quantile(target_pop, 0.025), # lower 95%CI
                tp_median = median(target_pop),
                tp_ub = quantile(target_pop, 0.975)) %>% # upper 95%CI
      mutate(ISO = country) %>%
      dplyr::select(ISO, year, tp_lb, tp_median, tp_ub)
    
    
    # cumulatively
    df_tp_cumu <- rc %>% 
      filter(is_target == 1) %>%
      group_by(year, run_id) %>%
      summarize(target_pop = sum(pop_model)) %>%
      group_by(run_id) %>%
      mutate(target_pop_cumu = cumsum(target_pop)) %>%
      dplyr::select(-target_pop) %>%
      group_by(year) %>%
      summarize(tp_cumu_lb = quantile(target_pop_cumu, 0.025), # lower 95%CI
                tp_cumu_median = median(target_pop_cumu),
                tp_cumu_ub = quantile(target_pop_cumu, 0.975)) %>% # upper 95%CI
      mutate(ISO = country) %>%
      dplyr::select(ISO, year, tp_cumu_lb, tp_cumu_median, tp_cumu_ub)
    
    # join by year and cumulative tables 
    df_tp <- left_join(df_tp_byyear, df_tp_cumu, by = c("ISO", "year"))
    
    # fill in no targeting years with 0s
    df_add <- data.frame("ISO" = country, "year" = vax_years[vax_years %in% df_tp$year == FALSE], 
                         "tp_lb" = 0, "tp_median" = 0, "tp_ub" = 0,
                         "tp_cumu_lb" = 0, "tp_cumu_median" = 0, "tp_cumu_ub" = 0) %>% as.data.frame()
    df_tp <- rbind(df_tp, df_add) %>% arrange(year)
    
    
    # add incidence columns
    df_case <- rc %>% 
      mutate(observed_case = pop_model * observed_incidence_rate) %>%
      group_by(year, run_id) %>%
      summarize(observed_case = sum(observed_case),
                true_case = sum(true_case)) %>%
      group_by(year) %>%
      summarize(observed_case_lb = quantile(observed_case, 0.025),
                observed_case_median = median(observed_case),
                observed_case_ub = quantile(observed_case, 0.975),
                true_case_lb = quantile(true_case, 0.025),
                true_case_median = median(true_case),
                true_case_ub = quantile(true_case, 0.975))
    
    pop_by_year <- rc %>%
      group_by(year, run_id) %>%
      summarize(pop=sum(pop_model)) %>%
      filter(run_id == 1) %>%
      dplyr::select(-run_id)
    
    df_inc <- df_case %>%
      left_join(pop_by_year, by = "year") %>%
      mutate(true_inc_rate_lb = true_case_lb / pop,
             true_inc_rate_median = true_case_median / pop,
             true_inc_rate_ub = true_case_ub / pop,
             obs_inc_rate_lb = observed_case_lb / pop,
             obs_inc_rate_median = observed_case_median / pop,
             obs_inc_rate_ub = observed_case_ub / pop) %>%
      dplyr::select(-c(true_case_lb, true_case_median, true_case_ub,
                       observed_case_lb, observed_case_median, observed_case_ub,
                       pop))
    
    # join df_tp with df_inc 
    df_result <- df_inc %>%
      left_join(df_tp, by = "year") %>%
      mutate(ISO = country) %>%
      dplyr::select(ISO, c(1:7, 9:14))
    
  }
  
  ### if want to get data for each admin
  if(sum_level == "admin"){
    
    if("NAME_2" %in% names(rc)){
      rc <- rc %>% dplyr::rename(admin_name = NAME_2)
      admin_level <- "admin2"
    }else{
      rc <- rc %>% dplyr::rename(admin_name = NAME_1)
      admin_level <- "admin1"
    }
    
    # target pop by year
    df_tp_byyear <- rc %>% 
      filter(is_target == 1) %>%
      group_by(admin_name, year) %>% 
      summarize(tp_lb = quantile(pop_model, 0.025), # lower 95%CI
                tp_median = median(pop_model),
                tp_ub = quantile(pop_model, 0.975)) %>% # upper 95%CI
      mutate(ISO = country) %>%
      dplyr::select(ISO, admin_name, year, tp_lb, tp_median, tp_ub) 

  
    # target pop cumulatively 
    df_tp_cumu <- rc %>% 
      filter(is_target == 1) %>%
      group_by(admin_name, year, run_id) %>%
      summarize(target_pop = sum(pop_model)) %>%
      group_by(run_id) %>%
      mutate(target_pop_cumu = cumsum(target_pop)) %>%
      dplyr::select(-target_pop) %>%
      group_by(admin_name, year) %>%
      summarize(tp_cumu_lb = quantile(target_pop_cumu, 0.025), # lower 95%CI
                tp_cumu_median = median(target_pop_cumu),
                tp_cumu_ub = quantile(target_pop_cumu, 0.975)) %>% # upper 95%CI
      mutate(ISO = country) %>%
      dplyr::select(ISO, admin_name, year, tp_cumu_lb, tp_cumu_median, tp_cumu_ub)
      
    
    # join by year and cumulative tables 
    df_tp <- left_join(df_tp_byyear, df_tp_cumu, by = c("ISO", "admin_name", "year"))
    
    # add back all the admin units/years that are not targeted, value = 0
    df_temp <- data.frame("admin_name" = rep(unique(rc$admin_name), times = length(vax_years)),
                          "year" = rep(vax_years, each=length(unique(rc$admin_name))))
    df_tp <- left_join(df_temp, df_tp, by = c("admin_name", "year")) %>%
      mutate(ISO = country) %>%
      dplyr::select(ISO, admin_name, year, tp_lb, tp_median, tp_ub, tp_cumu_lb, tp_cumu_median, tp_cumu_ub)
    
    # replace NA with 0
    df_tp[is.na(df_tp)] <-0
    
    
    ## add incidence columns
    df_inc <- rc %>%
      group_by(admin_name, year) %>%
      summarize(true_inc_rate_lb = quantile(true_incidence_rate, 0.025), 
                true_inc_rate_median = median(true_incidence_rate),
                true_inc_rate_ub = quantile(true_incidence_rate, 0.975),
                obs_inc_rate_lb = quantile(observed_incidence_rate, 0.025), 
                obs_inc_rate_median = median(observed_incidence_rate),
                obs_inc_rate_ub = quantile(observed_incidence_rate, 0.975))
    
    # join df_inc and df_tp
    df_result <- df_inc %>%
      left_join(df_tp, by = c("admin_name", "year")) %>%
      mutate(ISO = country) %>%
      arrange(year, admin_name) %>%
      dplyr::select(ISO, admin_name, year, true_inc_rate_lb, true_inc_rate_median, true_inc_rate_ub, obs_inc_rate_lb, obs_inc_rate_median, obs_inc_rate_ub, tp_lb, tp_median, tp_ub, tp_cumu_lb, tp_cumu_median, tp_cumu_ub)
    colnames(df_result)[2] <- admin_level
    
  }
  
  # get number of doses administered
  df_result <- df_result %>%
    mutate(doses_lb = vax_cov * tp_lb * 2, doses_median = vax_cov * tp_median * 2, doses_ub = vax_cov * tp_ub * 2,
           doses_cumu_lb = vax_cov * tp_cumu_lb * 2, doses_cumu_median = vax_cov * tp_cumu_median * 2, doses_cumu_ub = vax_cov * tp_cumu_ub * 2)

  return(df_result)
}



# test get_tp_doses_helper()
# a <- get_inc_tp_doses_helper(rc = rc_admin1, sum_level = "country", surveillance_scenario = "district-estimate")
# a <- get_inc_tp_doses_helper(rc = rc_admin1, sum_level = "admin", surveillance_scenario = "district-estimate")
# a <- get_inc_tp_doses_helper(rc = rc_admin2, sum_level = "country", surveillance_scenario = "district-estimate")
# a <- get_inc_tp_doses_helper(rc = rc_admin2, sum_level = "admin", surveillance_scenario = "district-estimate")
