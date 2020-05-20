library(data.table)
library(fst)
library(ggplot2)
source("src/functions.R")

#
w <- 11 #plot width (inches)
h <- 6.9 #plot height (inches)

#
benchmark <- read_fst(file.path("data", "processed", "constant_benchmark.fst"), as.data.table = TRUE)
model <- read_fst(file.path("data", "processed", "model_benchmark.fst"), as.data.table = TRUE)
benchmark <- benchmark[!is.na(target), .(state, date, target, days_left = as.integer(days_left), ci_upper, ci_lower, predicted_deaths, SCRPS)]
benchmark <- benchmark[days_left != 14]

# temp (something is wrong here...)
setnames(model, c("state", "date"), c("date", "state"))
model[, state := state + 14]
model[, days_left := as.integer(state - date)]
model <- model[, .(state, date, target, days_left, ci_upper, ci_lower, predicted_deaths, SCRPS)]
# end temp

# Pick three dates at random to plot
set.seed(1234)
example_dates <- sample(seq(as.Date("2020-04-15"), model[days_left == 13, max(date)], 1), 4)
plot_data <- rbindlist(list("Historical Avg." = benchmark[date %in% example_dates], "Capture-Retain" = model[date %in% example_dates]), idcol = "type")

plot <- ggplot(data = plot_data, aes(x = factor(days_left), color = type, group = type)) +
    geom_hline(aes(yintercept = target), color = "grey50") +
    geom_line(aes(y = predicted_deaths)) +
    geom_point(aes(y = predicted_deaths)) +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.3) +
    facet_wrap(~date) +
    set_default_theme() +
    labs(title = "Predicting deaths for a specific date, reported within 14 days.",
         subtitle = "",
         caption = "Source: FHM.",
         color = "Model",
         x = "Days of lag to predict",
         y = "Number of deaths")

ggsave(filename = file.path("output", "plots", "lag_prediction_by_date.pdf"),
       plot = plot, device = cairo_pdf, width = w, height = h)


plot_data <- rbindlist(
    list("Historical Avg." =
        benchmark[, .("Mean SCRPS" = mean(SCRPS, na.rm=T),
                      "Mean CI width" = mean(ci_upper-ci_lower, na.rm = TRUE),
                      "Share of CI covers true value" = mean((target <= ci_upper) * (target >= ci_lower), na.rm = TRUE)),
                  by = .(days_left)]#,
        # "Capture-Retain" =
        # model[, .(SCRPS    = mean(SCRPS, na.rm=T),
        #           CIwidth  = mean(ci_upper-ci_lower, na.rm = TRUE),
        #           coverage = mean((target <= ci_upper) * (target >= ci_lower), na.rm = TRUE)),
        #       by = .(days_left)],
        ), idcol = "type")

plot_data <- melt(plot_data, id.vars = c("type", "days_left"))

plot <- ggplot(data = plot_data, aes(x = factor(days_left), color = type, group = type)) +
    geom_line(aes(y = value)) +
    geom_point(aes(y = value)) +
    facet_wrap(~variable, scales = "free_y") +
    set_default_theme() +
    labs(title = "Model metrics",
         subtitle = "",
         caption = "Source: FHM.",
         color = "Model",
         x = "Days of lag to predict",
         y = "")

ggsave(filename = file.path("output", "plots", "model_metrics.pdf"),
       plot = plot, device = cairo_pdf, width = w, height = h)
