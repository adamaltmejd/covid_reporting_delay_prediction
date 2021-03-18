# Required packages
library(data.table)
library(lubridate)
library(fst)
library(readxl)
library(stringr)

##
# run file to create and reproduce the figures in the article
# make sure that correct updated FHM files exist at data/FHM/
# new files can be downloaded from <https://github.com/adamaltmejd/covid>
##
source("src/data_processing.R")
source("src/buildData.R")
source("src/constant_benchmark.R")
source("src/Beta_GP_lag_benchmark_build.R")
source("src/Beta_GP_lag_benchmark_refit.R")
source("src/plots.R")
