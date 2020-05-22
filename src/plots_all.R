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


plot_data <- rbindlist(list("Historical Avg." = benchmark, "Capture-Retain" = model), idcol = "type")
plot_data[, type := factor(type)]

dates <- seq(as.Date("2020-04-15"), model[days_left == 13, max(date)], 1)
plot_data_full <- plot_data
for(i in 1:length(dates)){
    plot_data <- plot_data_full[date %in% dates[i]]

colors <- c("#ECCBAE", "#046C9A", "gray50")
colors <- setNames(colors, c(levels(plot_data$type), "Reported"))

plot <- ggplot(data = plot_data,
               aes(x = days_left,
                   y = predicted_deaths,
                   color = type,
                   group = type)) +
    geom_hline(data = plot_data[!is.na(target), .(unique(target)), by = .(date)],
               aes(yintercept = V1), color = "grey50") +
    geom_line(data = reported_dead[date %in% dates[i]],
              aes(y = reported_dead, group = "Reported", color = "Reported"), linetype = "dashed") +
    geom_point(data = reported_dead[date %in% dates[i]],
               aes(y = reported_dead, group = "Reported", color = "Reported")) +
    # The actual models
    geom_line(position = position_dodge(width = 0.6)) +
    geom_point(data = plot_data[days_left != "0"], position = position_dodge(width = 0.6)) +
    geom_errorbar(data = plot_data[days_left != "0"], aes(ymin = ci_lower, ymax = ci_upper),
                  width = 0.7, position = position_dodge(width = 0.6)) +
    # Converging with a grey point on the last day
    geom_point(data = plot_data[days_left == "0" & type == "Historical Avg."], color = "grey50") +
    set_default_theme() +
    scale_color_manual(values = colors) +
    scale_y_continuous(minor_breaks = seq(0,200,10), breaks = seq(0,200,40), expand = expansion(add = c(0, 5))) +
    labs(title =format(dates[[i]],"%Y-%m-%d"),
         subtitle = "",
         caption = "",
         color = "Model",
         x = "Days of lag to predict",
         y = "Number of deaths")
ggsave(filename = file.path("output", "plots", paste("predicition",format(dates[[i]],"%Y-%m-%d"),'.pdf',sep="")),
       plot = plot, device = cairo_pdf, width = w, height = h/1.9)
}
