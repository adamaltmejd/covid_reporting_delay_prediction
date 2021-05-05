##
# building base data file from FHM data files
#
# reads files:
#  data/FHM/Folkhalsomyndigheten_Covid19_YYYY-MM-DD.xlsx
# (external files avialable at <https://github.com/adamaltmejd/covid>)
#
# generates the files:
# /data/processed/deaths_dt_SWE.fst
##

library(data.table)
library(fst)
library(stringr)
library(readxl)
source("src/data_processing_functions.R")

# FHM Deaths
files <- list.files(file.path("data", "FHM"), full.names = TRUE)
dts <- lapply(files, load_fhm_deaths)
deaths_dt <- rbindlist(dts)
setkey(deaths_dt, publication_date, date)
deaths_dt <- prepare_dt(deaths_dt)
write_fst(deaths_dt, file.path("data", "processed", "deaths_dt_SWE.fst"))

###
# Death stats UK
deaths_dt_UK <- fread(file.path("data", "uk", "uk.csv"))
setkey(deaths_dt_UK, publication_date, date)
deaths_dt_UK <- prepare_dt(deaths_dt_UK)
write_fst(deaths_dt_UK, file.path("data", "processed", "deaths_dt_UK.fst"))

###
# ICU stats
dts <- lapply(files, load_fhm_icu)
icu_dt <- rbindlist(dts)
setkey(icu_dt, publication_date, date)
icu_dt <- prepare_dt(icu_dt)
write_fst(icu_dt, file.path("data", "processed", "icu_dt_SWE.fst"))

