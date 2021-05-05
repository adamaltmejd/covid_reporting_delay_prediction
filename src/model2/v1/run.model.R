graphics.off()

##
# data input

##

run_model_default <- function(state) {

    data <- readRDS(file.path("data", "processed", "processed_data.rds"))

    model_parameters <- list(sim           = 400,
                             burnin        = ceiling(0.5*sim),
                             N.days.fixed  =  3,
                             quantile      = c(0.1,0.9))

    prior_list <- list(mu_beta    = c(0,0,0),
                       Sigma_beta = 1/2*diag(3),
                       a_sigma    = c(3,3),
                       b_sigma    = c(5/2,5/2))

    return(run_model(data, state, model_parameters, prior_list))
}

run_model <- function(data, state, model_parameters, prior_list) {
    source(file.path("src", "model2", "v1", "model.R"))
    require(numDeriv)
    require(invgamma)

    j <- grep(state, data$dates) # ! Is this correct?

    start_ = j - 31 # run the last 31 days

    report_cleaned <- report_clean(data$detected[start_:j,start_:j],data$dates[start_:j])
    new_cases <- newCases(report_cleaned)

    rownames(new_cases) <- as.character(data$dates[start_:j])
    colnames(new_cases) <- as.character(data$dates[start_:j])

    model(new_cases, model_parameters, prior_list)
}

fetch_plot_data <- function(model_results, deaths_dt, target_lag = 14) {
    source("src/util/util.scores.R")
    out <- as.data.table(model_results$Npost)

    out[, date := as.Date(dates)]
    out[, dates := NULL]

    target <- deaths_dt[days_lag == target_lag, .(date, target = N)]

    # TODO: need to fix deaths_dt so that it contains all days lag, before it was collapsed by work day

    # merge(out, target, by = "date", all.x = TRUE)

    # deaths_dt[date == "2020-10-28"]
    # deaths_dt[date == "2020-11-12"]

    # out[, SCRPS := SCRPS()]

    # SCRPS <- function(y, mu, sigma)

}


a <- run_model_default("2020-11-26")
b <- run_model_default("2020-06-20")
deaths_dt <- read_fst("data/processed/deaths_dt_SWE.fst", as.data.table = TRUE)


