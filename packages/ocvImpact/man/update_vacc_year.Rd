% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/surveillance_vacc_targeting.R
\name{update_vacc_year}
\alias{update_vacc_year}
\title{update_vacc_year}
\usage{
update_vacc_year(
  datapath,
  modelpath,
  country,
  rc_list,
  model_year,
  campaign_cov,
  threshold,
  surveillance_scenario,
  vac_start_year,
  vac_end_year,
  num_skip_years,
  rc_targeted,
  use_mean_ir,
  mean_ir_span,
  testing_sensitivity
)
}
\value{
rc_list
}
\description{
update the targets in the list if it's a vaccination campaign year
}
