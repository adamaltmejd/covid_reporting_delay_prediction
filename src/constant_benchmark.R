library(data.table)
library(fst)

generate_predictions <- function(deaths_dt) {
    # Calculate the average lag from the last 21 days for each report
    # So when calculating the average number of deaths added on day 3,
    # include day-3 reports from three weeks back from that date.
    DT <- deaths_dt[days_since_publication != 0 &
                   !is.na(days_since_publication) &
                   !is.na(date) &
                   date >= "2020-04-02",
                   .(date, publication_date, days_since_publication, n_diff)]

    # Create predictions for each publication date so we can track and evaluate
    # historical predictions
    report_dates <- seq(as.Date("2020-04-14"), deaths_dt[, max(publication_date)], 1)
    dts <- vector(mode = "list", length = length(report_dates))
    for (i in seq_along(report_dates)) {
        avg_delay <- DT[publication_date <= report_dates[i]]

        # For each delay day, calculate the mean additional deaths added
        # for the two weeks of reports preceding that date.
        # --> If we are calculating deaths added 7 days after, we look at
        #     deaths added 1 week ago back to deaths added 3 weeks ago.
        #     This way, we always take the mean of 2 weeks of reports (when available).
        avg_delay[, ref_date := max(date), by = days_since_publication]
        dts[[i]] <- avg_delay[date %between% list(ref_date - 14, ref_date),
                              .(avg_diff = mean(n_diff, na.rm = TRUE),
                                sd_diff = sd(n_diff, na.rm = TRUE)),
                              by = days_since_publication]
    }

    names(dts) <- report_dates
    avg_delay <- rbindlist(dts, idcol = "publication_date")
    avg_delay[, publication_date := as.Date(publication_date)]

    setkey(avg_delay, publication_date, days_since_publication)

    # To create actual predictions of totals per day
    # we need to add the averages to the reported data
    predictions <- avg_delay[deaths_dt[date >= "2020-04-02" & publication_date >= "2020-04-14"],
              on = .(publication_date,
                     days_since_publication > days_since_publication),
              by = .EACHI,
              .(date,
                sure_deaths = N,
                predicted_deaths = sum(avg_diff, na.rm = TRUE),
                predicted_deaths_SD = sqrt(sum(sd_diff^2, na.rm = TRUE)))] # assuming independently normal

    setnames(predictions, "publication_date", "prediction_date")
    predictions[, total := sure_deaths + predicted_deaths]

    # Assume no more deaths after 28 days just to have a cleaner data set
    predictions <- predictions[days_since_publication <= 28]
    predictions[, days_since_publication := NULL]

    setkey(predictions, prediction_date, date)

    return(predictions)
}

deaths_dt <- read_fst(file.path("data", "processed", "deaths_dt.fst"), as.data.table = TRUE)
DT <- generate_predictions(deaths_dt)
setkey(DT, prediction_date, date)
write_fst(DT, file.path("data", "processed", "constant_benchmark.fst"))
