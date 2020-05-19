library(data.table)
library(fst)
library(tidyr)
library(ggplot2)

generate_predictions <- function(deaths_dt) {
    # Calculate the average lag from the last 21 days for each report
    # So when calculating the average number of deaths added on day 3,
    # include day-3 reports from three weeks back from that date.
    DT <- deaths_dt[days_since_publication != 0 &
                   !is.na(days_since_publication) &
                   !is.na(date) &
                   date >= "2020-04-02",
                   .(date, publication_date, lag = days_since_publication, n_diff)]

    # To simplify, we assume no additional deaths are added after 28 days.
    DT <- DT[lag <= 28]

    # Create predictions for each publication date so we can track and evaluate
    # historical predictions
    report_dates <- seq(as.Date("2020-04-16"), deaths_dt[, max(publication_date)], 1)
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
    predictions <- avg_delay[deaths_dt[date >= "2020-04-02" & publication_date >= "2020-04-14"],
              on = .(state = publication_date), #lag > days_since_publication),
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

    # Don't include predictions states where we already know more than we are predicting
    predictions <- predictions[state - date < prediction_date - date + 1]

    # Drop unused columns
    predictions[, c("tmp.deaths_added", "tmp.deaths_added_sd") := NULL]

    return(predictions)
}

deaths_dt <- read_fst(file.path("data", "processed", "deaths_dt.fst"), as.data.table = TRUE)
DT <- generate_predictions(deaths_dt)

# For evaluation, try to predict the reported dead 14 days after the death date:
DT <- DT[prediction_date - date == 14]
DT[, target := reported_dead[state - date == 14], by = prediction_date]

write_fst(DT, file.path("data", "processed", "constant_benchmark.fst"))

plot_data <- DT[date == "2020-05-01"]
plot_data[, x := as.numeric(state - date)]
plot_data[, se := predicted_deaths_SD_cum / sqrt(14)]

g <- ggplot(data = plot_data, aes(x = x)) +
    geom_line(aes(y = target), color = "grey50") +
    geom_line(aes(y = predicted_deaths_cum)) +
    geom_point(aes(y = predicted_deaths_cum)) +
    geom_errorbar(aes(ymin = predicted_deaths_cum - 1.96 * se,
                      ymax = predicted_deaths_cum + 1.96 * se),
                  width = 0.1)

print(g)
