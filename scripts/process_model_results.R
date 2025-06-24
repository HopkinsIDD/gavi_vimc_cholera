# This script is used to process the final output from the VIMC core models
# Preamble ----------------------------------------------------------------

library(tidyverse)
library(ocvImpact)

# Load data ---------------------------------------------------------------
filepath <- 'COD_no-vaccination_zero_incid_trend_FALSE_outb_layer_FALSE_stoch.csv'
output <- read.csv(filepath)

# Adjust model outputs ----------------------------------------------------
adj_output <- output %>% 
  adjust_cases(., multiplier = 1/0.328) %>% 
  adjust_deaths(., multiplier  = 3.87) %>% 
  adjust_DALYs(., cases_multiplier = 1/0.328, deaths_multiplier = 3.87)
