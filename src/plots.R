##
# generate the figures for the article x
#
# read files:
# data/processed/deaths_dt.fst           (data_processing.R)
# data/processed/model_latest.fst        (Beta_GP_lag_benchmark_build.R)
# data/processed/constant_benchmark.fst  (constant_benchmark.R)
# data/processed/model_benchmark.fst.fst (Beta_GP_lag_benchmark_build.R)
#
# generates files:
# output/plots/SCRPS_over_weekdays.pdf
# output/plots/daily/prediction_YYYY-MM-DD.pdf
# output/plots/SCRPS_over_states.pdf
# output/plots/model_metrics.pdf
##

library(data.table)
library(fst)
library(forcats)
library(ggplot2)
library(hrbrthemes)
library(wesanderson)
source(file.path("src", "util","functions.R"))

#
w <- 11 # plot width (inches)
my_palette <- c("#d1ae90", "#046C9A", "#D69C4E", "#ABDDDE", "#000000") # wes_palette("Darjeeling2") # (replaced #ECCBAE)

# Plot 1 = Predictions and current stats
deaths_dt <- read_fst(file.path("data", "processed", "deaths_dt.fst"), as.data.table = TRUE)
model_predict <- read_fst(file.path("data", "processed", "model_latest.fst"), as.data.table = TRUE)
DT1 <- model_predict[date %between% c("2020-05-20", "2020-06-03")]
DT2 <- deaths_dt[!is.na(N) & !is.na(date) & publication_date == max(model_predict[, date])]
DT2[, avg := frollmean(N, 7, algo = "exact", align = "right")]
DT2 <- DT2[date > "2020-03-15"]

colors <- c("gray80", "#ECCBAE")
colors <- setNames(colors, c("Reported dead", "Model prediction"))

plot <- ggplot(data = DT1, aes(x = date, y = predicted_deaths)) +
    geom_bar(aes(fill = "Model prediction"), stat = "identity") +
    geom_bar(data = DT2,
             aes(y = N, fill = "Reported dead"), stat = "identity") +
    geom_line(data = DT2[!is.na(avg) & date <= max(model_predict[, date]) - 10],
              aes(y = avg, linetype = "7-day moving average"), color = "grey50") +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper),
                  width = 0.4, color = "#e0ac7e") +
    set_default_theme() +
    scale_fill_manual(values = colors) +
    scale_x_date(date_breaks = "3 day", date_labels = "%b %d", expand = expansion(add = 0.8)) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1.3, vjust = 1.1)) +
    scale_y_continuous(minor_breaks = seq(0,200,10), breaks = seq(0,200,40), expand = expansion(add = c(0, 5))) +
    labs(#title = paste0("Reported deaths as of ", deaths_dt[, max(publication_date)], " and model prediction"),
         #subtitle = "",
         #caption = "",
         fill = "",
         linetype = "",
         x = "Death date",
         y = "Number of deaths")

ggsave(filename = file.path("output", "paper", "plots", "latest_prediction.pdf"),
       plot = plot, device = cairo_pdf, width = w, height = w/1.9)


# Load data
benchmark <- read_fst(file.path("data", "processed", "constant_benchmark.fst"), as.data.table = TRUE)
model <- read_fst(file.path("data", "processed", "model_benchmark.fst"), as.data.table = TRUE)
benchmark <- benchmark[state >="2020-04-21"]

# Fix data tables so they look the same
reported_dead <- benchmark[, .(state, date, days_left = as.integer(days_left), reported_dead)]
benchmark <- benchmark[!is.na(target), .(state, date, target, days_left = as.integer(days_left), ci_upper, ci_lower, predicted_deaths, SCRPS)]
model[, days_left := as.integer(days_left)]
model[days_left == 0, `:=`(ci_upper = NA_real_, ci_lower = NA_real_, SCRPS = NA_real_)]

# Add a final 14-day prediction to model (just equal to truth)
# model <- rbindlist(list(model, data.table(state = model[, unique(date)] + 14, date = model[, unique(date)], days_left = 14)), use.names = TRUE, fill = TRUE)

# Days left as factor
reported_dead[, days_left := forcats::fct_rev(factor(days_left))]
benchmark[, days_left := forcats::fct_rev(factor(days_left))]
model[, days_left := forcats::fct_rev(factor(days_left))]

# Order correctly
setkey(reported_dead, date, state, days_left)
setkey(benchmark, date, state, days_left)
setkey(model, date, state, days_left)

# Add target as prediction for last day
model[is.na(target), predicted_deaths := model[!is.na(target), unique(target), by = date][, V1]]

#
## PLOT 1: Performance per day, PLOT1=4 random dates ##
#

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
        geom_point(data = DT[days_left != "0"], position = position_dodge(width = 0.6)) +
        geom_errorbar(data = DT[days_left != "0"], aes(ymin = ci_lower, ymax = ci_upper),
                    width = 0.7, position = position_dodge(width = 0.6)) +
        # Converging with a grey point on the last day
        geom_point(data = DT[days_left == "0" & type == "Historical Avg."], color = "grey50") +
        # Theming
        set_default_theme() +
        # scale_fill_manual(values = fill_colors, limits = label_order, drop = FALSE) +
        # scale_color_manual(values = my_palette) +
        scale_color_manual(values = colors) +
        scale_y_continuous(minor_breaks = seq(0,200,10), breaks = seq(0,200,40), expand = expansion(add = c(1, 5))) +
        labs(title = plot.title,
             #subtitle = "",
             #caption = "",
            color = "Model",
            x = "Days of lag to predict",
            y = "Number of deaths")

    if (DT[, uniqueN(date)] > 1) {
        # One plot for each day
        plot <- plot + facet_wrap(~date)
    }
    return(plot)
}

plot_data <- rbindlist(list("Historical Avg." = benchmark, "Capture-Retain" = model), idcol = "type")
plot_data[, type := factor(type)]

# Figure 1 - Pick 4 dates at random to plot
set.seed(1234)
example_dates <- sample(seq(as.Date("2020-04-21"), model[days_left == 13, max(date)], 1), 4)
plot <- day_plot(plot_data[date %in% example_dates],
                 reported_dead[date %in% example_dates],
                 "")
                 #"Predicting the number of deaths in a given day reported within 14 days.")
ggsave(filename = file.path("output", "paper", "plots", "lag_prediction_by_date.pdf"),
       plot = plot, device = cairo_pdf, width = w, height = w)

# For verification, plot each date as well
dates <- seq(as.Date("2020-04-26"), model[days_left == 13, max(date)], 1)

for (i in seq_along(dates)) {
    plot <- day_plot(plot_data[date == dates[i]],
                     reported_dead[date == dates[i]],
                     plot.title = dates[i])
    ggsave(filename = file.path("output", "paper", "plots", "daily", paste0("prediction_", dates[i], ".pdf")),
           plot = plot, device = cairo_pdf, width = w, height = w)
}
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
    scale_color_manual(values = my_palette) +
    labs(#title = "Model metrics",
         #subtitle = "",
         #caption = "",
         color = "Model",
         x = "Days of lag to predict",
         y = "")

ggsave(filename = file.path("output", "paper", "plots", "model_metrics.pdf"),
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
    scale_color_manual(values = my_palette) +
    labs(#title = "Model metrics",
         #subtitle = "",
         #caption = "",
         color = "Model",
         x = "Last date included in the model",
         y = "SCRPS")

ggsave(filename = file.path("output", "paper", "plots", "SCRPS_over_states.pdf"),
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
    scale_color_manual(values = my_palette) +
    labs(#title = "Model metrics",
         #subtitle = "",
         #caption = "",
         color = "Model",
         x = "Weekday",
         y = "")

ggsave(filename = file.path("output", "paper", "plots", "SCRPS_over_weekdays.pdf"),
       plot = plot, device = cairo_pdf, width = w, height = w/1.9)
