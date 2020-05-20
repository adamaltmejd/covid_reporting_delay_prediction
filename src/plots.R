library(data.table)
library(fst)
library(forcats)
library(ggplot2)
library(hrbrthemes)
library(wesanderson)
source("src/functions.R")

#
w <- 11 #plot width (inches)
h <- 6.9 #plot height (inches)

# Load data
benchmark <- read_fst(file.path("data", "processed", "constant_benchmark.fst"), as.data.table = TRUE)
model <- read_fst(file.path("data", "processed", "model_benchmark.fst"), as.data.table = TRUE)

# Fix data tables so they look the same
benchmark <- benchmark[!is.na(target), .(state, date, target, days_left = as.integer(days_left), ci_upper, ci_lower, predicted_deaths, SCRPS)]
model[, days_left := as.integer(days_left)]
# Add a final 14-day prediction to model (just equal to truth)
model <- rbindlist(list(model, data.table(state = model[, unique(date)] + 14, date = model[, unique(date)], days_left = 14)), use.names = TRUE, fill = TRUE)

# Days left as factor
benchmark[, days_left := forcats::fct_rev(factor(14 - days_left))]
model[, days_left := forcats::fct_rev(factor(14 - days_left))]

# Order correctly
setkey(benchmark, date, state, days_left)
setkey(model, date, state, days_left)

# Add target as prediction for last day
model[is.na(target), predicted_deaths := model[!is.na(target), unique(target), by = date][, V1]]

plot_data <- rbindlist(list("Historical Avg." = benchmark, "Capture-Retain" = model), idcol = "type")

# Pick three dates at random to plot
set.seed(1234)
example_dates <- sample(seq(as.Date("2020-04-15"), model[days_left == 13, max(date)], 1), 4)
plot_data <- plot_data[date %in% example_dates]

plot <- ggplot(data = plot_data, aes(x = factor(days_left), y = predicted_deaths, color = type, group = type)) +
    geom_hline(data = plot_data[!is.na(target), .(unique(target)), by = .(date)], aes(yintercept = V1), color = "grey50") +
    geom_line(alpha = 0.7) +
    geom_point(data = plot_data[days_left != "0"], alpha = 0.7) +
    geom_point(data = plot_data[days_left == "0" & type == "Historical Avg."], color = "grey50") +
    geom_errorbar(data = plot_data[days_left != "0"], aes(ymin = ci_lower, ymax = ci_upper), width = 0.4) +
    facet_wrap(~date) +
    set_default_theme() +
    scale_color_manual(values = wes_palette("Darjeeling2")) +
    labs(title = "Predicting the number of deaths for four dates reported within 14 days.",
         subtitle = "",
         caption = "Source: FHM.",
         color = "Model",
         x = "Days of lag to predict",
         y = "Number of deaths")

ggsave(filename = file.path("output", "plots", "lag_prediction_by_date.pdf"),
       plot = plot, device = cairo_pdf, width = w, height = h)


plot_data <- rbindlist(
    list("Historical Avg." =
         benchmark[, .(mean(SCRPS, na.rm=T),
                       mean(ci_upper-ci_lower, na.rm = TRUE),
                       mean((target <= ci_upper) * (target >= ci_lower), na.rm = TRUE)),
                  by = .(days_left)],
        "Capture-Retain" =
        model[, .(mean(SCRPS, na.rm=T),
                  mean(ci_upper-ci_lower, na.rm = TRUE),
                  mean((target <= ci_upper) * (target >= ci_lower), na.rm = TRUE)),
              by = .(days_left)]
        ), idcol = "type")

setnames(plot_data, c("V1", "V2", "V3"),
         c("Mean SCRPS", "Mean CI width (numer of deaths)", "Share of CI covers true value"))
plot_data <- plot_data[days_left != "0"]

plot_data <- melt(plot_data, id.vars = c("type", "days_left"))

plot <- ggplot(data = plot_data, aes(x = factor(days_left), color = type, group = type)) +
    geom_line(aes(y = value)) +
    geom_point(aes(y = value)) +
    facet_wrap(~variable, scales = "free_y") +
    set_default_theme() +
    scale_color_manual(values = wes_palette("Darjeeling2")) +
    labs(title = "Model metrics",
         subtitle = "",
         caption = "Source: FHM.",
         color = "Model",
         x = "Days of lag to predict",
         y = "")

ggsave(filename = file.path("output", "plots", "model_metrics.pdf"),
       plot = plot, device = cairo_pdf, width = w, height = h)
