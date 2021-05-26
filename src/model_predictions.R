library(fst)
library(data.table)

max.days.to.report <- 30

# Prepare datasets
source(file.path("src", "buildData.R"))
buildData("sweden", file.path("data", "processed", "processed_data_sweden.rds"))
buildData("uk", file.path("data", "processed", "processed_data_uk.rds"))

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

source(file.path("src","util","util_swe.r"))
swe_data <- readRDS(file.path("data", "processed", "processed_data_sweden.rds"))

#remove data before  2020-07-01
start_date <- as.Date("2020-07-01")
index <- swe_data$dates_report >= start_date
swe_data$detected             = swe_data$detected[index,index]
swe_data$report               = swe_data$report[index,index]
swe_data$dates                = swe_data$dates[index]
swe_data$dates_report         = swe_data$dates_report[index]
swe_data$dates_not_reported   = swe_data$dates_not_reported[index]

target_swe <- data.frame(reported = swe_data$detected[row(swe_data$detected)+max.days.to.report==col(swe_data$detected)])
target_swe$dates <- swe_data$dates_report[1:length(target_swe$reported)]
#remove obs above max.days to report
swe_data$detected[row(swe_data$detected)+max.days.to.report<col(swe_data$detected)]=NA

dts_swe <- lapply(
    as.Date(swe_data$dates_report[swe_data$dates_report >= start_date + 90]),
    FUN = function(x, ...) swe.prediction(report.dates = x, ...),
    max.days.to.report = max.days.to.report,
    result = swe_data,
    target = target_swe
)
<<<<<<< HEAD
dts_smooth_swe <- lapply(dts, gp.smooth, max.days.to.report = max.days.to.report)
write_fst(rbindlist(dts), file.path("data", "model_predictions_full_SWE.fst"))
write_fst(rbindlist(dts_smooth_swe), file.path("data", "model_predictions_full_SWE_smooth.fst"))
=======
write_fst(rbindlist(dts_swe), file.path("data", "processed", "model_predictions_full_SWE.fst"))

dts_smooth_swe <- lapply(dts_swe, gp.smooth, max.days.to.report = max.days.to.report)
write_fst(rbindlist(dts_smooth_swe), file.path("data", "processed", "model_predictions_full_smooth_SWE.fst"))
>>>>>>> 368df4d84405f5780bb5370388558059698984d4

################################
# UK Prediction over all dates #
################################

source(file.path("src", "util", "util_uk.r"))
uk_data <- readRDS(file.path("data", "processed", "processed_data_uk.rds"))
#zero.report <- uk_data$dates %in% as.Date(c("2021-01-26", "2021-01-28", "2021-03-01"))
#uk_data$report[, zero.report] <- 0

target_uk <- data.frame(reported = uk_data$detected[row(uk_data$detected)+max.days.to.report==col(uk_data$detected)])
target_uk$dates <- uk_data$dates_report[1:length(target_uk$reported)]
uk_data$detected[row(uk_data$detected)+max.days.to.report<col(uk_data$detected)]=NA

dts_uk <- lapply(
    as.Date(uk_data$dates_report[uk_data$dates_report >= "2020-10-10"]),
    FUN = function(x, ...) uk.prediction(report.dates = x, ...),
    max.days.to.report = max.days.to.report,
    result = uk_data,
    target = target_uk
)
<<<<<<< HEAD
dts_smooth <- lapply(dts_uk, gp.smooth, max.days.to.report = max.days.to.report)
write_fst(rbindlist(dts_smooth), file.path("data", "model_predictions_full_UK_smooth.fst"))
write_fst(rbindlist(dts_uk), file.path("data", "model_predictions_full_UK.fst"))
#res <- unlist(lapply(dts, function(x){x[date==state+4,(ci_lower<= target)* (ci_upper>= target)]}))
=======
write_fst(rbindlist(dts_uk), file.path("data", "processed", "model_predictions_full_UK.fst"))

dts_smooth_uk <- lapply(dts_uk, gp.smooth, max.days.to.report = max.days.to.report)
write_fst(rbindlist(dts_smooth_uk), file.path("data", "processed", "model_predictions_full_smooth_UK.fst"))
>>>>>>> 368df4d84405f5780bb5370388558059698984d4
