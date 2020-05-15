library(data.table)
library(fst)

DT <- read_fst(file.path("data", "processed", "deaths_dt.fst"), as.data.table = TRUE)

# Calculate the average lag from the last 21 days for each report
# So when calculating the average number of deaths added on day 3,
# include day-3 reports from three weeks back from that date.

# Calculate the predictions for the last 3 weeks to get a track record back in time.
days <- 21
dts <- vector(mode = "list", length = days)
for (i in 1:days) {
    avg_delay <- DT[days_since_publication != 0 &
                    !is.na(days_since_publication) &
                    !is.na(date) &
                    date >= "2020-04-02"]
    avg_delay[, ref_date := max(date) + 1 - i, by = days_since_publication]
    dts[[i]] <- avg_delay[
        date %between% list(ref_date - 21, ref_date),
        .(avg_diff = mean(n_diff, na.rm = TRUE)),
        by = days_since_publication]
}

names(dts) <- DT[days_since_publication == 1, max(date)] - 20:0

avg_delay <- rbindlist(dts, idcol = "report_date")
avg_delay[, report_date := as.Date(report_date)]

write_fst(death_dt, file.path("data", "processed", "constant_benchmark.fst"))
