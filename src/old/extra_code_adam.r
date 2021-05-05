deaths_dt


deaths_dt <- read_fst(file.path("data", "processed", "deaths_dt_SWE.fst"), as.data.table = TRUE)

# Include only up to 14 days
n_days <- 14
DT <- deaths_dt[days_lag <= n_days]
DT <- DT[date > "2020-04-02"]
DT <- DT[date <= max(publication_date) - n_days]
# DT[, N_14d := N[days_lag == 14], by = date]

# Ignore negative reports
DT[n_diff < 0, n_diff := 0]
DT <- DT[!is.na(workdays_lag)]

lag_means <- function(lags, N, groups = seq(0, 1, 0.25)) {
    means <- data.table(cut(1:sum(N), quantile(1:sum(N), probs = groups), include = TRUE), rep(lags, N))[, .(l=max(V2)), by = V1][, l]
    out <- data.table(t(means))
    setnames(out, paste0("p", gsub("\\.", "", groups[2:length(groups)])))

    return(out)
}

plot_data <- melt(DT[, lag_means(as.integer(workdays_lag), n_diff, groups = seq(0, 1, 0.25)), by = date], id.vars = "date")

g <- ggplot(data = plot_data, aes(x = date, y = value, group = variable, color = variable)) + geom_line()
g

model[date == "2020-05-20", ]

g <- ggplot(data = deaths_dt[date > "2020-04-02", sum(n_diff), by = .(date, workdays_lag <= 2)], aes(x = date, y = V1, color = workdays_lag, group = workdays_lag)) + geom_point() + geom_smooth()
g

deaths_dt[date > "2020-04-02", sum(n_diff), by = .(date, workdays_lag <= 1)][date <= "2020-04-15"]
deaths_dt[date == "2020-04-12"]
