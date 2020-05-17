library(data.table)
library(fst)

calculate_average_lags <- function(deaths_dt) {
    # Calculate the average lag from the last 21 days for each report
    # So when calculating the average number of deaths added on day 3,
    # include day-3 reports from three weeks back from that date.

    # Calculate the predictions for the last 3 weeks to get a track record back in time.
    days <- 21
    dts <- vector(mode = "list", length = days)
    for (i in 1:days) {
        avg_delay <- deaths_dt[days_since_publication != 0 &
                               !is.na(days_since_publication) &
                               !is.na(date) &
                               date >= "2020-04-02"]
        avg_delay[, ref_date := max(date) + 1 - i, by = days_since_publication]
        dts[[i]] <- avg_delay[
            date %between% list(ref_date - 21, ref_date),
            .(avg_diff = mean(n_diff, na.rm = TRUE), sd_diff = sd(n_diff, na.rm = TRUE)),
            by = days_since_publication]
    }

    names(dts) <- deaths_dt[days_since_publication == 1, max(date)] - 20:0

    avg_delay <- rbindlist(dts, idcol = "publication_date")
    avg_delay[, publication_date := as.Date(publication_date)]
    setkey(avg_delay, publication_date, days_since_publication)

    return(avg_delay)
}

generate_predictions <- function(deaths_dt, avg_delay) {
    # To create actual predictions of totals per day
    # we need to add the averages to the reported data
    avg_delay[deaths_dt,
              on = .(publication_date,
                     days_since_publication > days_since_publication),
              by = .EACHI,
              .(date, N,
                pred_N = N + sum(avg_diff, na.rm = TRUE),
                pred_SD = sqrt(sum(sd_diff^2, na.rm = TRUE)))] # assuming independently normal
}

deaths_dt <- read_fst(file.path("data", "processed", "deaths_dt.fst"), as.data.table = TRUE)
avg_delay <- calculate_average_lags(deaths_dt)

DT <- generate_predictions(deaths_dt, avg_delay)
setkey(DT, publication_date, date)

write_fst(DT, file.path("data", "processed", "constant_benchmark.fst"))
