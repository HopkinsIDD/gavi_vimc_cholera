#### This script adds documentation to package data ####

#' @collate data_documentation.R


#' surveillance_confirmation_rate_parameters
#'
#' Parameters for beta distribution that fits the confirmation rate distribution, Used for the surveillance project
#' 
#'
#' @format ## `surveillance_confirmation_rate_parameters`
#' A data frame with rows and columns:
#' \describe{
#'   \item{shape1}{parameter 1}
#'   \item{shape2}{parameter 2}
#'   ...
#' }

#' default_country_incidence
#'
#' Country name, country code, and incidence rate for 47 VIMC countries for which we do not have modeled raster estimates. 
#' These country-level estimates are projected to a grid and used to run the VIMC model for these countries. Last updated in 2020
#'
#' @format ## `default_country_incidence`
#' A data frame with rows and columns:
#' \describe{
#'   \item{country}{country name}
#'   \item{country_code}{country code}
#'   \item{singular_estimate}{singular estimate}
#'   \item{is_num_case}{either 1 or NA}
#'   \item{year_list}{year range for estimate}
#'   \item{data_source}{data source}
#'   \item{num_case}{number of cases}
#'   \item{incid_rate_100k}{incidence rate per 100k people}
#'   ...
#' }


#' who_annual_reports
#'
#' Reported data on country level suspected cases and deaths, used to calibrate incidence rate trends. source: WHO Annual Cholera Reports, last updated in 2020
#'
#' @format ## `who_annual_reports`
#' A data frame with ? rows and ? columns:
#' \describe{
#'   \item{ISO3}{country ISO3 code}
#'   \item{Country}{country name}
#'   \item{Continent}{continent}
#'   \item{Region1}{region in continent}
#'   \item{TL}{time left}
#'   \item{TR}{time right}
#'   \item{year}{year}
#'   \item{sCh}{number of suspected cases}
#'   \item{deaths}{number of deaths}
#'   ...
#' }


#' Annual case fatality ratios by country. Calculated from case and death data in WHO Annual Cholera Reports. Last updated in 2019.
#'
#' @format ## `who_cfrs`
#' A data frame with 85 rows and 7 columns:
#' \describe{
#'   \item{year}{year}
#'   \item{cntry_code}{country code}
#'   \item{cases}{number of cases}
#'   \item{deaths}{number of deaths}
#'   \item{cfr}{case fatality ratio}
#'   \item{country_name}{country name}
#'   \item{source}{source}
#'   ...
#' }

#' Modeled one dose vaccine effectiveness decay function fit from only observational studies, used in the 2023 touchstone, source: Xu and Azman 2024
#'
#' @format ## `ve_log_1dose_obs`
#' A :
#' \describe{
#'   ...
#' }

#' Modeled two dose vaccine effectiveness decay function fit from only observational studies, used in the 2023 touchstone, source: Xu and Azman 2024
#'
#' @format ## `ve_log_2dose_obs`
#' A :
#' \describe{
#'   ...
#' }

#' Disability weights used to calculate DALYs and YLLs in 2019, 2021, and 2023 touchstones, source: IHME Global Burden of Disease, 2019.
#'
#' @format ## `disability_weights`
#' A data frame with 2120 rows and 4 columns:
#' \describe{
#'   \item{Sequela}{sequela}
#'   \item{Health state name}{health state}
#'   \item{Health state lay description}{description of health state}
#'   \item{Disability Weight}{disability weight}
#'   ...
#' }
#' @source IHME Global Burden of Disease, 2019

#' List of default countries that could be modeled in the VIMC core model, as called in scripts/set_all_parameters.R
#' The default value is currently overwritten by another list of countries hard-coded in set_all_parameters
#'
#' @format ## `default_countries`
#' A data frame with 45 rows and 1 columns:
#' \describe{
#'   \item{country}{country code}
#'   ...
#' }

#' ## Custom coverage table used for the 2024 DRC vaccine impact case study using the format used by VIMC
#' This represents an example file for scenarios with "custom" coverage settings (when custom$use_montagu_coverage: FALSE, custom$coverage_filename)
#'
#' @format ## `drc_custom_coverage_2024_2026`
#' A data frame with 6 rows and 15 columns:
#' \describe{
#'   \item{scenario}{vaccination scenario}
#'   \item{set_name}{one dose, two dose vaccination set name}
#'   \item{vaccine}{either OCV1 or OCV2}
#'   \item{gavi_support}{how much of this campaign is supported by GAVI}
#'   \item{activity_type}{type of vaccination, campaign for cholera}
#'   \item{country_code}{country code}
#'   \item{country}{country name}
#'   \item{year}{year campaign will take place}
#'   \item{age_first}{age of youngest vaccinee}
#'   \item{age_last}{age of oldest vaccinee}
#'   \item{age_range_verbatim}{age range}
#'   \item{target}{target population size}
#'   \item{coverage}{proportion of the target population covered by vaccine}
#'   \item{gender}{which gender is vaccinated}
#'   \item{proportion_risk}{proportion at risk}
#'   ...
#' }




