runname <- "201910gavi-5"
scenarios <- c("campaign-default", "no-vaccination")
# countries <- c("COD", "ETH", "KEN", "SOM", "SSD", "HTI")
countries <-    c("AGO", "BDI", "BEN", "BFA", "CAF", "CIV", "CMR", "COD", "COG", "DZA", "ETH", "GHA", "GIN", "GNB", "KEN", "LBR", "MDG", "MLI",
                  "MOZ", "MRT", "MWI", "NAM", "NER", "NGA", "RWA", "SEN", "SLE", "SOM", "SSD", "TCD", "TGO", "TZA", "UGA", "ZAF", "ZMB", "ZWE",
                  "AFG", "HTI", "IRN", "IRQ", "NPL", "PAK", "PHL", "THA", "YEM", "IND", "BGD")
clean_outputs <- TRUE
targeting_strategy <- "affected_pop"
clean_incid <- FALSE
num_samples <- 30
num_skip_years <- 3