library(fst)
library(data.table)

###
# Sweden predictions over all dates

###
# UK Prediction over all dates
source(file.path("src", "util", "util_uk.r"))
uk_data <- readRDS(file.path("data", "processed", "processed_data_uk.rds"))

zero.report <- uk_data$dates %in% as.Date(c("2021-01-26", "2021-01-28", "2021-03-01"))
max.days.to.report <- 30
uk_data$report[, zero.report] <- 0

dts <- lapply(
    as.Date(uk_data$dates_report[uk_data$dates_report >= "2020-10-10"]),
    FUN = function(x, ...) uk.prediction(report.dates = x, ...),
    max.days.to.report = max.days.to.report,
    result = uk_data
)
dts_smooth <- lapply(dts, gp.smooth, max.days.to.report = max.days.to.report)

write_fst(rbindlist(dts_smooth), file.path("data", "model_predictions_full_uk.fst"))
