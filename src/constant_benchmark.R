##
# builds the constant benchmark model
#
# load files:
# data/processed/deaths_dt_SWE.fst (data_processing.R)
#
# generates the files:
# data/processed/constant_benchmark.fst
##

library(data.table)
library(fst)
library(tidyr)
library(ggplot2)
source(file.path("src","util","util.scores.R"))

alpha <- 0.1

generate_predictions <- function(deaths_dt) {
    # Calculate the average lag from the last 21 days for each report
    # So when calculating the average number of deaths added on day 3,
    # include day-3 reports from three weeks back from that date.
    DT <- deaths_dt[days_lag != 0 &
                    !is.na(days_lag) &
                    !is.na(date) &
                    date >= "2020-04-02",
                    .(date, publication_date, lag = days_lag, n_diff)]

    # To simplify, we assume no additional deaths are added after 30 days.
    DT <- DT[lag <= 30]

    # Create predictions for each publication date so we can track and evaluate
    # historical predictions
    report_dates <- seq(DT[, min(publication_date)] + 30, DT[, max(publication_date)], 1)
    dts <- vector(mode = "list", length = length(report_dates))
    for (i in seq_along(report_dates)) {
        tmp.delay <- DT[publication_date <= report_dates[i]]

        # For each delay day, calculate the mean additional deaths added
        # for the two weeks of reports preceding that date.
        # --> If we are calculating deaths added 7 days after, we look at
        #     deaths added 1 week ago back to deaths added 3 weeks ago.
        #     This way, we always take the mean of 2 weeks of reports (when available).
        tmp.delay[, ref_date := max(date), by = lag]
        tmp.delay <- tmp.delay[date %between% list(ref_date - 13, ref_date),
                               .(avg_diff = mean(n_diff, na.rm = TRUE),
                                 sd_diff = sd(n_diff, na.rm = TRUE),
                                 n_obs = sum(!is.na(n_diff))),
                               by = lag]

        dts[[i]] <- tmp.delay
    }

    names(dts) <- report_dates
    avg_delay <- rbindlist(dts, idcol = "state")
    avg_delay[, state := as.Date(state)]
    setkey(avg_delay, state, lag)

    # Then create a table of predictions, where we calculate the expected number
    # of deaths reported for a death_date at a given report_date.
    # So for example, with the _state_ of knowledge at date 2020-05-16,
    # we want predict the number of deaths that happened on 2020-05-10 and have
    # been reported until 2020-05-11, 12, etc.
    predictions <- avg_delay[deaths_dt[date >= "2020-04-02" & publication_date >= "2020-04-30"],
              on = .(state = publication_date), #lag > days_lag),
              by = .EACHI,
              .(date,
                prediction_date = date + 1:.N,
                reported_dead = N,
                #sd_n_obs = n_obs
                tmp.deaths_added = tidyr::replace_na(avg_diff, 0),
                tmp.deaths_added_sd = tidyr::replace_na(sd_diff, 0))]

    setkey(predictions, state, date, prediction_date)

    # Don't predict anything for dates we know the truth about
    predictions[prediction_date <= state, `:=`(tmp.deaths_added = 0,
                                               tmp.deaths_added_sd = 0)]
    predictions[, `:=`(predicted_deaths_cum = reported_dead + cumsum(tmp.deaths_added),
                       predicted_deaths_SD_cum = sqrt(cumsum(tmp.deaths_added_sd^2))), # assuming independently normal
                by = .(state, date)]
    #predictions[predicted_deaths_SD_cum == 0, predicted_deaths_SD_cum := NA_real_]

    # Don't include predictions states where we already know more than we are predicting
    predictions <- predictions[state - date < prediction_date - date + 1]

    # Drop unused columns
    predictions[, c("tmp.deaths_added", "tmp.deaths_added_sd") := NULL]

    return(predictions)
}

##########
# SWEDEN #
##########

deaths_dt_SWE <- read_fst(file.path("data", "processed", "deaths_dt_SWE.fst"), as.data.table = TRUE)
DT <- generate_predictions(deaths_dt_SWE)

# For evaluation, try to predict the reported dead 30 days after the death date:
DT <- DT[prediction_date - date == 30]
DT[, target := reported_dead[state - date == 30], by = prediction_date]

out <- DT[, .(state,
              date,
              predicted_deaths = predicted_deaths_cum,
              predicted_deaths_SD = predicted_deaths_SD_cum,
              reported_dead,
              target,
              days_left = 30 - as.integer(state - date),
              ci_upper = predicted_deaths_cum + qnorm(1-alpha/2) * predicted_deaths_SD_cum,
              ci_lower = predicted_deaths_cum + qnorm(alpha/2) * predicted_deaths_SD_cum,
              CRPS = CRPS(target, predicted_deaths_cum, predicted_deaths_SD_cum))]

out[is.nan(CRPS), CRPS := 0]
write_fst(out, file.path("data", "processed", "constant_model_predictions_full_SWE.fst"))

######
# UK #
######

deaths_dt_UK <- read_fst(file.path("data", "processed", "deaths_dt_UK.fst"), as.data.table = TRUE)
DT <- generate_predictions(deaths_dt_UK)

# For evaluation, try to predict the reported dead 30 days after the death date:
DT <- DT[prediction_date - date == 30]
DT[, target := reported_dead[state - date == 30], by = prediction_date]

out <- DT[, .(state,
              date,
              predicted_deaths = predicted_deaths_cum,
              predicted_deaths_SD = predicted_deaths_SD_cum,
              reported_dead,
              target,
              days_left = 30 - as.integer(state - date),
              ci_upper = predicted_deaths_cum + qnorm(1-alpha/2) * predicted_deaths_SD_cum,
              ci_lower = predicted_deaths_cum + qnorm(alpha/2) * predicted_deaths_SD_cum,
              CRPS = CRPS(target, predicted_deaths_cum, predicted_deaths_SD_cum))]
out[is.nan(CRPS), CRPS := 0]

write_fst(out, file.path("data", "processed", "constant_model_predictions_full_UK.fst"))
