##
# run file to create and reproduce the figures in the article
# make sure that correct FHM files exists at data/FHM/
# they are aviable from <https://github.com/adamaltmejd/covid>
##
source("src/data_processing.R")
source("src/buildData.R")
source("src/constant_benchmark.R")
source("src/Beta_GP_lag_benchmark_build.R")
source("src/Beta_GP_lag_benchmark_refit.R")
source("src/plots.R")
