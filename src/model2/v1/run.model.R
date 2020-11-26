graphics.off()
source(file.path("src","model2","v1","model.R"))
library(numDeriv)
library(invgamma)
##
# data input
## 
sim <- 40000
data <- readRDS(file.path("data", "processed", "processed_data.rds"))
j <- dim(data$detected)[1]
start_ = j - 31 # run the last 31 days
data_cut <- data
data_cut$detected           <- data_cut$detected[start_:j,start_:j]
data_cut$dates              <- data_cut$dates[start_:j]
data_cut$dates_report       <-  data_cut$dates_report[start_:j]
data_cut$dates_not_reported <- data_cut$dates_not_reported[start_:j]

model_parameters <- list(sim           = sim,
                         burnin        = ceiling(0.5*sim),
                         N.days.fixed  =  1,
                         quantile      = c(0.1,0.9))

prior_list <- list(mu_beta    = c(0,0,0),
                   Sigma_beta = 1/2*diag(3),
                   a_sigma    = c(3,3),
                   b_sigma    = c(5/2,5/2))

result <- model(data_cut, model_parameters, prior_list)