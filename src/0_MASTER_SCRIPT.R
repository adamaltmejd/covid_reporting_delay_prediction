##
# run file to create and reproduce the figures in the article
# make sure that correct updated FHM files exist at data/FHM/
# new files can be downloaded from <https://github.com/adamaltmejd/covid>
##
cat(paste0("\n[", Sys.time(), "]", "running file: src/data_processing.R\n"))
source("src/data_processing.R")
cat(paste0("[", Sys.time(), "]", "running file: src/constant_benchmark.R\n"))
source("src/constant_benchmark.R")
cat(paste0("[", Sys.time(), "]", "running file: src/model_predictions.R\n"))
source("src/model_predictions.R")
cat(paste0("[", Sys.time(), "]", "running file: src/plots.R\n"))
source("src/plots.R")
cat(paste0("[", Sys.time(), "]", "DONE.\n"))
