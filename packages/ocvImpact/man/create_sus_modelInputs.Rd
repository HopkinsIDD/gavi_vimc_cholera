% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_sus_modelInputs.R
\name{create_sus_modelInputs}
\alias{create_sus_modelInputs}
\title{create_sus_modelInputs}
\usage{
create_sus_modelInputs(
  datapath,
  modelpath,
  country,
  scenario,
  rawoutpath,
  vacc_alloc,
  cache,
  clean,
  trunc_year = 10
)
}
\arguments{
\item{datapath}{path to input data}

\item{modelpath}{path to montagu files}

\item{country}{country code}

\item{scenario}{Unique string that identifies the coverage scenario name}

\item{rawoutpath}{path to raw model output files}

\item{vacc_alloc}{object returned from [`allocate_vaccine()`]}

\item{cache}{montagu cache}

\item{clean}{logical that indicates whether existing sus_files should be deleted}

\item{trunc_year}{number of years after which VE = 0 (default = 10)}
}
\value{
track_prop_immune dataframe
}
\description{
Create raster for proportion of population susceptible in each year, after accounting for vaccination
}
