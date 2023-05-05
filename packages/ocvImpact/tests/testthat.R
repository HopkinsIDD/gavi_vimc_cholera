library(testthat)
oxygen2::roxygenise("packages/ocvImpact")
install.packages("packages/ocvImpact", type = "source", repos = NULL)
library('ocvImpact', character.only = T)

test_local(path = "packages/ocvImpact")
