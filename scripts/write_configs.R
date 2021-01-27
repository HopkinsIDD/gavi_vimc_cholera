# roxygen2::roxygenise("packages/ocvImpact")
# install.packages("packages/ocvImpact", type = "source", repos = NULL)

source("scripts/set_all_parameters.R")

cpathname <- file.path("configs", runname)
dir.create(cpathname, showWarnings = FALSE)

for(scn in scenarios){

  scnpathname <- file.path(cpathname, scn)
  dir.create(scnpathname, showWarnings = FALSE)
  pars <- tidyr::expand_grid(runname = runname, scenario = scn, country = countries, targeting = targeting_strategy, nsamples = num_samples, nskipyear = num_skip_years, clean = clean_outputs, redrawIncid = clean_incid)

  lapply(1:nrow(pars), function(i){
    par_row <- pars[i,]
    ocvImpact::prepare_config(par_row, scnpathname)
  })

}



