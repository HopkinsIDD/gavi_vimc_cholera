% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_montagu.R
\name{import_centralburden_template}
\alias{import_centralburden_template}
\title{import_centralburden_template}
\usage{
import_centralburden_template(modelpath, country, cache, redownload = TRUE)
}
\arguments{
\item{modelpath}{Path to montagu files}

\item{country}{country code}

\item{cache}{montagu cache}

\item{redownload}{whether to redownload the file}
}
\value{
dataframe for central burden template for one country
}
\description{
Imports a CSV template file for the central burden estimates static columns indicating the burden estimates that need to be generated for each scenario
}
