% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/incidence_rate_trend_multiplier.R
\name{incidence_rate_trend_multiplier}
\alias{incidence_rate_trend_multiplier}
\title{incidence_rate_trend_multiplier}
\usage{
incidence_rate_trend_multiplier(
  datapath,
  modelpath,
  country,
  datathreshold = 5,
  use_country_incid_trend
)
}
\arguments{
\item{datapath}{path to input data}

\item{modelpath}{path to montagu files}

\item{country}{country code}

\item{datathreshold}{the threshold for the minimum number of years of data in the WHO dataset to give a valid incidence rate estimate (the default is 5)}

\item{use_country_incid_trend}{whether we wanna use the country-level who data for a certain country. It overrides others.}
}
\value{
a function that can take in years and give a incidence rate prediction for a given country
}
\description{
Create and return a function that can take in the year to give a prediction of the incidence rate in that country
}
