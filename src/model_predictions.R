library(fst)
library(data.table)

# Creates .fst files containing data tables with model predictions for all different states.
# Each file should have the following columns:
# * state: date at which the data was recorded
# * date: the date
# * days_left: state - date
# * predicted_deaths: the prediciton of the model
# * ci_lower: lower CI of prediction
# * ci_upper: upper CI of prediction
# * CRPS: accuracy score
# * target: the "true" number of deaths within 30 days (that the model was trained to predict)

#####################################
# Sweden predictions over all dates #
#####################################

# Prepare Swedish data
source(file.path("src", "buildData.R"))
buildData("sweden", file.path("data", "processed", "processed_data_sweden.rds"))
# Load prepared data
swe_data <- readRDS(file.path("data", "processed", "processed_data_sweden.rds"))

# Build Sweden model
# TODO
# zero.report <- swe_data$dates %in% as.Date(c("2021-01-26", "2021-01-28", "2021-03-01"))
# max.days.to.report <- 30
# swe_data$report[, zero.report] <- 0

#dt <- swe.prediction(report.dates = as.Date(swe_data$dates_report[swe_data$dates_report >= "2020-10-10"])[10], max.days.to.report = max.days.to.report, result = swe_data)
#dt_smooth <- gp.smooth(dt, max.days.to.report = max.days.to.report)
# dts <- lapply(
#     as.Date(swe_data$dates_report[swe_data$dates_report >= "2020-10-10"]),
#     FUN = function(x, ...) swe.prediction(report.dates = x, ...),
#     max.days.to.report = max.days.to.report,
#     result = swe_data
# )
# dts_smooth <- lapply(dts, gp.smooth, max.days.to.report = max.days.to.report)
# write_fst(rbindlist(dts_smooth), file.path("data", "model_predictions_full_SWE.fst"))

################################
# UK Prediction over all dates #
################################

# Prepare UK data
source(file.path("src", "buildData.R"))
source(file.path("src", "util", "util_uk.r"))
buildData("uk", file.path("data", "processed", "processed_data_uk.rds"))
# Load prepared data
uk_data <- readRDS(file.path("data", "processed", "processed_data_uk.rds"))

zero.report <- uk_data$dates %in% as.Date(c("2021-01-26", "2021-01-28", "2021-03-01"))
max.days.to.report <- 30
uk_data$report[, zero.report] <- 0


target <- data.frame(reported = result$detected[row(result$detected)+max.days.to.report==col(result$detected)])
target$dates <- result$dates_report[1:length(target$reported)]
#remove obs  above max.days to report
result$detected[row(result$detected)+max.days.to.report<col(result$detected)]=NA

#dt <- uk.prediction(report.dates = as.Date(uk_data$dates_report[uk_data$dates_report >= "2020-10-10"])[10], max.days.to.report = max.days.to.report, result = uk_data)
#dt_smooth <- gp.smooth(dt, max.days.to.report = max.days.to.report)

dts <- lapply(
    as.Date(uk_data$dates_report[uk_data$dates_report >= "2020-10-10"]),
    FUN = function(x, ...) uk.prediction(report.dates = x, ...),
    max.days.to.report = max.days.to.report,
    result = uk_data,
    target = target
)
dts_smooth <- lapply(dts, gp.smooth, max.days.to.report = max.days.to.report)
write_fst(rbindlist(dts_smooth), file.path("data", "model_predictions_full_UK.fst"))
