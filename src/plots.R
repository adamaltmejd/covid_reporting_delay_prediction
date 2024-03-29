##
# generate the figures for the article
#
# read files:
# data/processed/deaths_dt_SWE.fst       (data_processing.R)
# data/processed/model_latest.fst        (Beta_GP_lag_benchmark_build.R)
# data/processed/constant_benchmark.fst  (constant_benchmark.R)
# data/processed/model_benchmark.fst     (Beta_GP_lag_benchmark_build.R)
#
# generates files:
# output/plots/CRPS_over_weekdays.pdf
# output/plots/daily/prediction_YYYY-MM-DD.pdf
# output/plots/CRPS_over_states.pdf
# output/plots/model_metrics.pdf
##

library(data.table)
library(fst)
library(forcats)
library(ggplot2)
library(hrbrthemes)
library(wesanderson)
#library(extrafont)
library(Cairo)
library(cowplot)
set_null_device("cairo")
source(file.path("src", "util", "functions.R"))

#
my_palette <- c("#d1ae90", "#046C9A", "#D69C4E", "#ABDDDE", "#000000")

##########
# PLOT 2: LATEST PREDICTIONS

latest_prediction <- function(deaths_dt, model_predict) {
    # Predict over the last 15 days, plot last two months
    plot_state <- model_predict[, max(state)]
    model_predict <- model_predict[state == plot_state & date %between% c(plot_state - 15, plot_state)]

    deaths_dt[, avg := frollmean(N, 7, algo = "exact", align = "center"), by = publication_date]
    deaths_dt <- deaths_dt[publication_date == plot_state & date %between% c(plot_state - 60, plot_state)]

    colors <- c("gray80", "#ECCBAE")
    colors <- setNames(colors, c("Reported dead", "Model prediction"))

    ggplot(data = model_predict, aes(x = date, y = predicted_deaths)) +
        geom_bar(aes(fill = "Model prediction"), stat = "identity") +
        geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper),
                      width = 0.4, color = "#e0ac7e") +
        geom_bar(data = deaths_dt,
                aes(y = N, fill = "Reported dead"), stat = "identity") +
        geom_line(data = deaths_dt[date <= plot_state - 10],
                  aes(y = avg, linetype = "7-day moving average"), color = "grey50") +
        set_default_theme() +
        scale_fill_manual(values = colors) +
        scale_x_date(date_breaks = "1 week", date_labels = "%b %d", expand = expansion(add = 0.8)) +
        theme(axis.text.x = element_text(angle = 35, hjust = 1.3, vjust = 1.1)) +
        scale_y_continuous(minor_breaks = seq(0,200,10), breaks = seq(0,200,10), expand = expansion(add = c(0, 5))) +
        labs(fill = "",
             linetype = "",
             x = "Death date",
             y = "Number of deaths")

}
deaths_dt_SWE <- read_fst(file.path("data", "processed", "deaths_dt_SWE.fst"), as.data.table = TRUE)
model_predict_SWE <- read_fst(file.path("data", "processed", "model_predictions_full_smooth_SWE.fst"), as.data.table = TRUE)
plot_SWE  <- latest_prediction(deaths_dt_SWE, model_predict_SWE)

deaths_dt_UK <- read_fst(file.path("data", "processed", "deaths_dt_UK.fst"), as.data.table = TRUE)
model_predict_UK <- read_fst(file.path("data", "processed", "model_predictions_full_smooth_UK.fst"), as.data.table = TRUE)
plot_UK <- latest_prediction(deaths_dt_UK, model_predict_UK)

p <- plot_grid(plot_grid(
    plot_SWE + guides(fill = "none", linetype = "none") + theme(plot.margin = unit(c(25, 5, 5, 5), "pt")),
    plot_UK + guides(fill = "none", linetype = "none") + ylab(NULL) + theme(plot.margin = unit(c(25, 5, 5, 5), "pt")),
    labels = c("Sweden", "United Kingdom"), label_fontfamily = "EB Garamond",
    hjust = -0.5, align = "hv",
    nrow = 1, ncol = 2),
    get_legend(plot_SWE),
    nrow = 2, ncol = 1, rel_heights = c(0.9, 0.1)
)
save_plot(filename = file.path("output", "paper", "plots", "latest_prediction.pdf"),
          plot = p, ncol = 2, nrow = 1, base_height = 4, device = cairo_pdf)

##########
# PREPARE DATA

# Load data
min_date <- as.Date("2020-10-10")
model_SWE <- read_fst(file.path("data", "processed", "model_predictions_full_smooth_SWE.fst"), as.data.table = TRUE)
max_state <- model_SWE[, max(state)]
model_SWE <- model_SWE[date >= min_date & state <= date + 30 & state <= max_state - 30]
model_SWE[, days_left := 30 - as.integer(state - date)]
benchmark_SWE <- read_fst(file.path("data", "processed", "constant_model_predictions_full_SWE.fst"), as.data.table = TRUE)
benchmark_SWE <- benchmark_SWE[date >= min_date & state <= date + 30 & state <= max_state - 30]
benchmark_SWE[days_left == 0, `:=`(ci_upper = NA_real_, ci_lower = NA_real_, CRPS = NA_real_)]
model_SWE[days_left == 0, `:=`(ci_upper = NA_real_, ci_lower = NA_real_, CRPS = NA_real_)]

model_UK <- read_fst(file.path("data", "processed", "model_predictions_full_smooth_UK.fst"), as.data.table = TRUE)
max_state <- model_UK[, max(state)]
model_UK <- model_UK[date >= min_date & state <= date + 30 & state <= max_state - 30]
model_UK[, days_left := 30 - as.integer(state - date)]
benchmark_UK <- read_fst(file.path("data", "processed", "constant_model_predictions_full_UK.fst"), as.data.table = TRUE)
benchmark_UK <- benchmark_UK[date >= min_date & state <= date + 30 & state <= max_state - 30]
benchmark_UK[days_left == 0, `:=`(ci_upper = NA_real_, ci_lower = NA_real_, CRPS = NA_real_)]
model_UK[days_left == 0, `:=`(ci_upper = NA_real_, ci_lower = NA_real_, CRPS = NA_real_)]

reported_dead_SWE <- benchmark_SWE[, .(state, date, days_left, reported_dead)]
reported_dead_UK <- benchmark_UK[, .(state, date, days_left, reported_dead)]

# Select same columns
model_SWE <- model_SWE[, .(state, date, days_left, target, predicted_deaths, ci_lower, ci_upper, CRPS)]
benchmark_SWE <- benchmark_SWE[, .(state, date, days_left, target, predicted_deaths, ci_lower, ci_upper, CRPS)]
model_UK <- model_UK[, .(state, date, days_left, target, predicted_deaths, ci_lower, ci_upper, CRPS)]
benchmark_UK <- benchmark_UK[, .(state, date, days_left, target, predicted_deaths, ci_lower, ci_upper, CRPS)]

# Order correctly
setkey(reported_dead_SWE, date, state, days_left)
setkey(reported_dead_UK, date, state, days_left)
setkey(benchmark_SWE, date, state, days_left)
setkey(benchmark_UK, date, state, days_left)
setkey(model_SWE, date, state, days_left)
setkey(model_UK, date, state, days_left)

##########
# PLOT 3: Performance per day for 3 random dates

# Plot a specific day
day_plot <- function(DT, reported, plot.title) {
    colors <- c("#d1ae90", "#046C9A", "gray50")
    colors <- setNames(colors, c(levels(DT$type), "Reported"))

    if (any(DT[!is.na(target), uniqueN(target) != 1, by = date][, V1])) {
        warning("Multiple unique target values for same day: ", plot.title)
    }
    hline <- DT[!is.na(target), .(unique(target)), by = .(date)]

    plot <- ggplot(data = DT,
                aes(x = days_left,
                    y = predicted_deaths,
                    color = type,
                    group = type)) +
        geom_hline(data = hline, aes(yintercept = V1), color = "grey50") +
        geom_line(data = reported,
                  aes(y = reported_dead, group = "Reported", color = "Reported"), linetype = "dashed") +
        geom_point(data = reported,
                   aes(y = reported_dead, group = "Reported", color = "Reported")) +
        # The actual models
        geom_line(aes(linetype = type), position = position_dodge(width = 0.6)) +
        geom_point(data = DT[days_left != 0], position = position_dodge(width = 0.6)) +
        geom_errorbar(data = DT[days_left != 0], aes(ymin = ci_lower, ymax = ci_upper),
                    width = 0.7, position = position_dodge(width = 0.6)) +
        # Converging with a grey point on the last day
        geom_point(data = DT[days_left == 0 & type == "Benchmark model"], color = "grey50") +
        # Theming
        set_default_theme() + guides(linetype = "none") +
        scale_color_manual(values = colors) +
        scale_y_continuous(breaks = scales::extended_breaks(), expand = expansion(add = c(1, 5))) +
        scale_x_reverse(breaks = scales::extended_breaks(), expand = expansion(add = c(1, 1))) +
        labs(title = plot.title,
             color = "Model",
             x = "Days of lag to predict", # lag prediction distance (days)
             y = "Number of deaths")

    if (DT[, uniqueN(date)] > 1) {
        # One plot for each day
        plot <- plot + facet_wrap(~date, scales = "free")
    }
    return(plot)
}

# PLOT DATA
plot_data_SWE <- rbindlist(list("Benchmark model" = benchmark_SWE, "Prediction model" = model_SWE), idcol = "type")
plot_data_SWE[, type := factor(type)]
plot_data_UK <- rbindlist(list("Benchmark model" = benchmark_UK, "Prediction model" = model_UK), idcol = "type")
plot_data_UK[, type := factor(type)]

# Figure 1 - Pick 4 dates at random to plot
set.seed(1234)
example_dates <- sample(seq(model_UK[, min(date)], model_SWE[, max(date)], 1), 3)
plot_SWE <- day_plot(plot_data_SWE[date %in% example_dates],
                     reported_dead_SWE[date %in% example_dates],
                     "")
# ggsave2(filename = file.path("output", "paper", "plots", "lag_prediction_by_date_SWE.pdf"),
#        plot = plot, device = cairo_pdf, width = 7, height = 7)

set.seed(1234)
example_dates <- sample(seq(model_UK[, min(date)], model_UK[, max(date)], 1), 3)
plot_UK <- day_plot(plot_data_UK[date %in% example_dates],
                    reported_dead_UK[date %in% example_dates],
                    "")
# ggsave2(filename = file.path("output", "paper", "plots", "lag_prediction_by_date_UK.pdf"),
#        plot = plot, device = cairo_pdf, width = 7, height = 7)

# Figure 1 - 3 dates per country
p <- plot_grid(plot_grid(
    plot_SWE + guides(fill = "none", linetype = "none", color = "none") + xlab(NULL),
    plot_UK + guides(fill = "none", linetype = "none", color = "none"),
    labels = c("Sweden", "United Kingdom"), label_fontfamily = "EB Garamond",
    hjust = 0, align = "hv",
    nrow = 2, ncol = 1),
    get_legend(plot_SWE),
    nrow = 2, ncol = 1, rel_heights = c(0.9, 0.1)
)
save_plot(filename = file.path("output", "paper", "plots", "lag_prediction_by_date.pdf"),
          plot = p, ncol = 1, nrow = 2, base_height = 5, device = cairo_pdf)

# For verification, plot each date as well
dates <- seq(model_SWE[, min(date)], model_SWE[, max(date) - 30], 1)
for (i in seq_along(dates)) {
    plot <- day_plot(plot_data_SWE[date == dates[i]],
                     reported_dead_SWE[date == dates[i]],
                     plot.title = dates[i])
    ggsave2(filename = file.path("output", "paper", "plots", "daily", "SWE", paste0("prediction_", dates[i], ".pdf")),
            plot = plot, device = cairo_pdf, width = 7, height = 7)
}

dates <- seq(model_UK[, min(date)], model_UK[, max(date) - 30], 1)
for (i in seq_along(dates)) {
    plot <- day_plot(plot_data_UK[date == dates[i]],
                     reported_dead_UK[date == dates[i]],
                     plot.title = dates[i])
    ggsave2(filename = file.path("output", "paper", "plots", "daily", "UK", paste0("prediction_", dates[i], ".pdf")),
            plot = plot, device = cairo_pdf, width = 7, height = 7)
}

###########
# PLOT 4: MODEL METRICS

plot_data_SWE <- rbindlist(
    list("Benchmark model" =
         benchmark_SWE[, .(mean(CRPS, na.rm=T),
                       mean(ci_upper-ci_lower, na.rm = TRUE),
                       mean((target <= ci_upper) * (target >= ci_lower), na.rm = TRUE)),
                  by = .(days_left)],
        "Prediction model" =
        model_SWE[, .(mean(CRPS, na.rm=T),
                  mean(ci_upper-ci_lower, na.rm = TRUE),
                  mean((target <= ci_upper) * (target >= ci_lower), na.rm = TRUE)),
              by = .(days_left)]
        ), idcol = "type")

setnames(plot_data_SWE, c("V1", "V2", "V3"),
         c("Mean CRPS", "Mean CI width (numer of deaths)", "Share of CI covers true value"))
plot_data_SWE <- plot_data_SWE[days_left != 0]

plot_data_SWE <- melt(plot_data_SWE, id.vars = c("type", "days_left"))

plot_SWE <- ggplot(data = plot_data_SWE, aes(x = days_left, color = type, group = type)) +
    geom_line(aes(y = value)) +
    geom_point(aes(y = value)) +
    facet_wrap(~variable, scales = "free_y") +
    set_default_theme() +
    scale_color_manual(values = my_palette) +
    scale_x_reverse(breaks = scales::extended_breaks(), expand = expansion(add = c(1, 1))) +
    scale_y_continuous(labels = function(x) scales::number(x, style_negative = "minus")) +
    labs(color = "Model",
         x = "Days of lag to predict",
         y = "")

# ggsave2(filename = file.path("output", "paper", "plots", "model_metrics_SWE.pdf"),
#        plot = plot, device = cairo_pdf, width = 5, height = 5/1.9)

plot_data_UK <- rbindlist(
    list("Benchmark model" =
         benchmark_UK[, .(mean(CRPS, na.rm=T),
                       mean(ci_upper-ci_lower, na.rm = TRUE),
                       mean((target <= ci_upper) * (target >= ci_lower), na.rm = TRUE)),
                  by = .(days_left)],
        "Prediction model" =
        model_UK[, .(mean(CRPS, na.rm=T),
                  mean(ci_upper-ci_lower, na.rm = TRUE),
                  mean((target <= ci_upper) * (target >= ci_lower), na.rm = TRUE)),
              by = .(days_left)]
        ), idcol = "type")

setnames(plot_data_UK, c("V1", "V2", "V3"),
         c("Mean CRPS", "Mean CI width (numer of deaths)", "Share of CI covers true value"))
plot_data_UK <- plot_data_UK[days_left != 0]

plot_data_UK <- melt(plot_data_UK, id.vars = c("type", "days_left"))

plot_UK <- ggplot(data = plot_data_UK, aes(x = days_left, color = type, group = type)) +
    geom_line(aes(y = value)) +
    geom_point(aes(y = value)) +
    facet_wrap(~variable, scales = "free_y") +
    set_default_theme() +
    scale_color_manual(values = my_palette) +
    scale_x_reverse(breaks = scales::extended_breaks(), expand = expansion(add = c(1, 1))) +
    scale_y_continuous(labels = function(x) scales::number(x, style_negative = "minus")) +
    labs(color = "Model",
         x = "Days of lag to predict",
         y = "")

# ggsave2(filename = file.path("output", "paper", "plots", "model_metrics_UK.pdf"),
#        plot = plot, device = cairo_pdf, width = 5, height = 5/1.9)

# Figure 1 - 3 dates per country
p <- plot_grid(plot_grid(
    plot_SWE + guides(fill = "none", linetype = "none", color = "none") + xlab(NULL),
    plot_UK + guides(fill = "none", linetype = "none", color = "none") + theme(plot.margin = unit(c(15, 5, 5, 5), "pt")),
    labels = c("Sweden", "United Kingdom"), label_fontfamily = "EB Garamond",
    hjust = 0, align = "hv",
    nrow = 2, ncol = 1),
    get_legend(plot_SWE),
    nrow = 2, ncol = 1, rel_heights = c(0.9, 0.1)
)
save_plot(filename = file.path("output", "paper", "plots", "model_metrics.pdf"),
          plot = p, ncol = 1, nrow = 2, base_height = 5, device = cairo_pdf)

###########
#  PLOT A1: Performance over time (as more training data becomes availiable)

# Only include dates for which we have all 30 dates
states <- fintersect(benchmark_SWE[!is.na(CRPS), .N, state][N == 30, .(state)],
                     model_SWE[!is.na(CRPS), .N, state][N == 30, .(state)])[ , state]

plot_data <- rbindlist(list(
    "Benchmark model" = benchmark_SWE[state %in% states & days_left != 0, mean(CRPS), by = state],
    "Prediction model" = model_SWE[state %in% states & days_left != 0, mean(CRPS), by = state]
    ), idcol = "type")

plot_SWE <- ggplot(data = plot_data, aes(x = state, y = V1, color = type, group = type)) +
    geom_line() + geom_point() +
    scale_y_continuous(limits = c(-80, 0), labels = function(x) scales::number(x, style_negative = "minus")) +
    set_default_theme() +
    scale_color_manual(values = my_palette) +
    labs(color = "Model",
         x = "",
         y = "Mean CRPS")

# ggsave2(filename = file.path("output", "paper", "plots", "CRPS_over_states_SWE.pdf"),
#        plot = plot, device = cairo_pdf, width = 5, height = 5/1.9)

states <- fintersect(benchmark_UK[!is.na(CRPS), .N, state][N == 30, .(state)],
                     model_UK[!is.na(CRPS), .N, state][N == 30, .(state)])[ , state]

plot_data <- rbindlist(list(
    "Benchmark model" = benchmark_UK[state %in% states & days_left != 0, mean(CRPS), by = state],
    "Prediction model" = model_UK[state %in% states & days_left != 0, mean(CRPS), by = state]
    ), idcol = "type")

plot_UK <- ggplot(data = plot_data, aes(x = state, y = V1, color = type, group = type)) +
    geom_line() + geom_point() +
    scale_y_continuous(limits = c(-80, 0), labels = function(x) scales::number(x, style_negative = "minus")) +
    set_default_theme() +
    scale_color_manual(values = my_palette) +
    labs(color = "Model",
         x = "Last date included in the model",
         y = "CRPS")

# ggsave2(filename = file.path("output", "paper", "plots", "CRPS_over_states_UK.pdf"),
#        plot = plot, device = cairo_pdf, width = 5, height = 5/1.9)

p <- plot_grid(plot_grid(
    plot_SWE + guides(fill = "none", linetype = "none", color = "none") + theme(plot.margin = unit(c(25, 5, 5, 5), "pt")),
    plot_UK + guides(fill = "none", linetype = "none", color = "none") + ylab(NULL) + theme(plot.margin = unit(c(25, 5, 5, 5), "pt")),
    labels = c("Sweden", "United Kingdom"), label_fontfamily = "EB Garamond",
    hjust = 0, align = "hv",
    nrow = 1, ncol = 2),
    get_legend(plot_SWE),
    nrow = 2, ncol = 1, rel_heights = c(0.9, 0.1)
)
save_plot(filename = file.path("output", "paper", "plots", "CRPS_over_states.pdf"),
          plot = p, ncol = 2, nrow = 1, base_width = 5, device = cairo_pdf)


###########
#  PLOT A2:  Performance by day-of-week

# Base results on same dates
dates <- fintersect(benchmark_SWE[!is.na(CRPS), .(date, state)], model_SWE[!is.na(CRPS), .(date, state)])
plot_data <- rbindlist(list("Benchmark model" = benchmark_SWE[dates],
                            "Prediction model" = model_SWE[dates]), idcol = "type")
plot_data[, dayofweek := factor(weekdays(date),
                                levels = c("Monday", "Tuesday", "Thursday",
                                           "Wednesday", "Friday", "Saturday",
                                           "Sunday"))]
plot_data <- plot_data[days_left != 0, mean(CRPS), by = .(dayofweek, type)]
plot_SWE <- ggplot(data = plot_data, aes(x = dayofweek, y = V1, color = type, group = type)) +
    geom_line() + geom_point() +
    scale_y_continuous(limits = c(-25, 0), labels = function(x) scales::number(x, style_negative = "minus")) +
    set_default_theme() +
    scale_color_manual(values = my_palette) +
    labs(color = "Model",
         x = "Weekday",
         y = "Mean CRPS")
# ggsave2(filename = file.path("output", "paper", "plots", "CRPS_over_weekdays_SWE.pdf"),
#        plot = plot, device = cairo_pdf, width = 5, height = 5/1.9)

dates <- fintersect(benchmark_UK[!is.na(CRPS), .(date, state)], model_UK[!is.na(CRPS), .(date, state)])
plot_data <- rbindlist(list("Benchmark model" = benchmark_UK[dates],
                            "Prediction model" = model_UK[dates]), idcol = "type")
plot_data[, dayofweek := factor(weekdays(date),
                                levels = c("Monday", "Tuesday", "Thursday",
                                           "Wednesday", "Friday", "Saturday",
                                           "Sunday"))]
plot_data <- plot_data[days_left != 0, mean(CRPS), by = .(dayofweek, type)]
plot_UK <- ggplot(data = plot_data, aes(x = dayofweek, y = V1, color = type, group = type)) +
    geom_line() + geom_point() +
    scale_y_continuous(limits = c(-25, 0), labels = function(x) scales::number(x, style_negative = "minus")) +
    set_default_theme() +
    scale_color_manual(values = my_palette) +
    labs(color = "Model",
         x = "Weekday",
         y = "Mean CRPS")
# ggsave2(filename = file.path("output", "paper", "plots", "CRPS_over_weekdays_UK.pdf"),
#        plot = plot, device = cairo_pdf, width = 5, height = 5/1.9)

p <- plot_grid(plot_grid(
    plot_SWE + guides(fill = "none", linetype = "none", color = "none") + theme(plot.margin = unit(c(25, 5, 5, 5), "pt")),
    plot_UK + guides(fill = "none", linetype = "none", color = "none") + ylab(NULL) + theme(plot.margin = unit(c(25, 5, 5, 5), "pt")),
    labels = c("Sweden", "United Kingdom"), label_fontfamily = "EB Garamond",
    hjust = 0, align = "hv",
    nrow = 1, ncol = 2),
    get_legend(plot_SWE),
    nrow = 2, ncol = 1, rel_heights = c(0.9, 0.1)
)
save_plot(filename = file.path("output", "paper", "plots", "CRPS_over_weekdays.pdf"),
          plot = p, ncol = 2, nrow = 1, base_width = 5, device = cairo_pdf)

