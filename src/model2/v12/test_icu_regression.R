###
# test icu_regression
#
###

source(file.path("src","model2","v12","icu_regression.R"))

deaths <- readRDS(file.path("data", "processed", "processed_data.rds"))
icu    <- readRDS(file.path("data", "processed", "processed_data_icu.rds"))
cov_ <- icu_covariates(deaths, icu)
