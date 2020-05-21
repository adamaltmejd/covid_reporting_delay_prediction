library(data.table)
library(fst)
library(forcats)
library(ggplot2)
library(hrbrthemes)
library(wesanderson)
source("src/functions.R")

#
w <- 11 # plot width (inches)

# Load data
benchmark <- read_fst(file.path("data", "processed", "constant_benchmark.fst"), as.data.table = TRUE)
model <- read_fst(file.path("data", "processed", "model_benchmark.fst"), as.data.table = TRUE)

# Fix data tables so they look the same
reported_dead <- benchmark[, .(state, date, days_left = as.integer(days_left), reported_dead)]
benchmark <- benchmark[!is.na(target), .(state, date, target, days_left = as.integer(days_left), ci_upper, ci_lower, predicted_deaths, SCRPS)]
model[, days_left := as.integer(days_left)]
# Add a final 14-day prediction to model (just equal to truth)
model <- rbindlist(list(model, data.table(state = model[, unique(date)] + 14, date = model[, unique(date)], days_left = 14)), use.names = TRUE, fill = TRUE)

# Days left as factor
reported_dead[, days_left := forcats::fct_rev(factor(14 - days_left))]
benchmark[, days_left := forcats::fct_rev(factor(14 - days_left))]
model[, days_left := forcats::fct_rev(factor(14 - days_left))]

# Order correctly
setkey(reported_dead, date, state, days_left)
setkey(benchmark, date, state, days_left)
setkey(model, date, state, days_left)

# Add target as prediction for last day
model[is.na(target), predicted_deaths := model[!is.na(target), unique(target), by = date][, V1]]

#
## PLOT 1: Performance on 4 random dates ##
#
plot_data <- rbindlist(list("Historical Avg." = benchmark, "Capture-Retain" = model), idcol = "type")
plot_data[, type := factor(type)]
# Pick three dates at random to plot
set.seed(1234)
example_dates <- sample(seq(as.Date("2020-04-15"), model[days_left == 13, max(date)], 1), 4)
plot_data <- plot_data[date %in% example_dates]

colors <- c("#ECCBAE", "#046C9A", "gray50")
colors <- setNames(colors, c(levels(plot_data$type), "Reported"))

plot <- ggplot(data = plot_data,
               aes(x = days_left,
                   y = predicted_deaths,
                   color = type,
                   group = type)) +
    geom_hline(data = plot_data[!is.na(target), .(unique(target)), by = .(date)],
               aes(yintercept = V1), color = "grey50") +
    geom_line(data = reported_dead[date %in% example_dates],
              aes(y = reported_dead, group = "Reported", color = "Reported"), linetype = "dashed") +
    geom_point(data = reported_dead[date %in% example_dates],
               aes(y = reported_dead, group = "Reported", color = "Reported")) +
    # The actual models
    geom_line(position = position_dodge(width = 0.6)) +
    geom_point(data = plot_data[days_left != "0"], position = position_dodge(width = 0.6)) +
    geom_errorbar(data = plot_data[days_left != "0"], aes(ymin = ci_lower, ymax = ci_upper),
                  width = 0.7, position = position_dodge(width = 0.6)) +
    # Converging with a grey point on the last day
    geom_point(data = plot_data[days_left == "0" & type == "Historical Avg."], color = "grey50") +
    # One plot for each day
    facet_wrap(~date) +
    # Theming
    set_default_theme() +
    # scale_fill_manual(values = fill_colors, limits = label_order, drop = FALSE) +
    # scale_color_manual(values = wes_palette("Darjeeling2")) +
    scale_color_manual(values = colors) +
    scale_y_continuous(minor_breaks = seq(0,200,10), breaks = seq(0,200,40), expand = expansion(add = c(0, 5))) +
    labs(title = "Predicting the number of deaths in a given day reported within 14 days.",
         subtitle = "",
         caption = "",
         color = "Model",
         x = "Days of lag to predict",
         y = "Number of deaths")

ggsave(filename = file.path("output", "plots", "lag_prediction_by_date.pdf"),
       plot = plot, device = cairo_pdf, width = w, height = w)

#
## PLOT 2: Statistics ##
#

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
         caption = "",
         color = "Model",
         x = "Days of lag to predict",
         y = "")

ggsave(filename = file.path("output", "plots", "model_metrics.pdf"),
       plot = plot, device = cairo_pdf, width = w, height = w/1.9)

#
## PLOT 3: Performance over time (as more training data becomes availiable) ##
#
# Only include dates for which we have all 14 dates
states <- fintersect(benchmark[!is.na(SCRPS), .N, state][N == 14, .(state)],
                     model[!is.na(SCRPS), .N, state][N == 14, .(state)])[ , state]

plot_data <- rbindlist(list(
    "Historical Avg." = benchmark[state %in% states & days_left != "0", mean(SCRPS), by = state],
    "Capture-Retain" = model[state %in% states & days_left != "0", mean(SCRPS), by = state]
    ), idcol = "type")

plot <- ggplot(data = plot_data, aes(x = state, y = V1, color = type, group = type)) +
    geom_line() + geom_point() +
    set_default_theme() +
    scale_color_manual(values = wes_palette("Darjeeling2")) +
    labs(title = "Model metrics",
         subtitle = "",
         caption = "",
         color = "Model",
         x = "Last date included in the model",
         y = "SCRPS")

ggsave(filename = file.path("output", "plots", "SCRPS_over_states.pdf"),
       plot = plot, device = cairo_pdf, width = w, height = w/1.9)

#
## PLOT 4: Performance by day-of-week ##
#
# Base results on same dates
dates <- fintersect(benchmark[!is.na(SCRPS), .(date, state)], model[!is.na(SCRPS), .(date, state)])
plot_data <- rbindlist(list("Historical Avg." = benchmark[dates],
                            "Capture-Retain" = model[dates]), idcol = "type")

plot_data[, dayofweek := factor(weekdays(date),
                                levels = c("Monday", "Tuesday", "Thursday",
                                           "Wednesday", "Friday", "Saturday",
                                           "Sunday"))]

plot_data <- plot_data[days_left != "0", mean(SCRPS), by = .(dayofweek, type)]

plot <- ggplot(data = plot_data, aes(x = dayofweek, y = V1, color = type, group = type)) +
    geom_line() + geom_point() +
    # facet_wrap(~days_left) +
    set_default_theme() +
    scale_color_manual(values = wes_palette("Darjeeling2")) +
    labs(title = "Model metrics",
         subtitle = "",
         caption = "",
         color = "Model",
         x = "Last date included in the model",
         y = "")

ggsave(filename = file.path("output", "plots", "SCRPS_over_weekdays.pdf"),
       plot = plot, device = cairo_pdf, width = w, height = h/1.9)
