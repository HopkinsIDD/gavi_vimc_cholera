# roxygen2::roxygenise("packages/ocvImpact")
# install.packages("packages/ocvImpact", type = "source", repos = NULL)

source("scripts/set_all_parameters.R")

pars <- tidyr::expand_grid(runname = runname, scenario = scenarios, country = countries, targeting = targeting_strategy, nsamples = num_samples, nskipyear = num_skip_years)
cpathname <- file.path("configs", runname)
dir.create(cpathname, showWarnings = FALSE)

lapply(1:nrow(pars), function(i){
  par_row <- pars[i,]
  ocvImpact::prepare_config(par_row, cpathname)
})

