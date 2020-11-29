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

report_cleaned <- report_clean(data$detected[start_:j,start_:j],data$dates[start_:j])
new_cases <- newCases(report_cleaned)
rownames(new_cases) <- as.character(data$dates[start_:j])
colnames(new_cases) <- as.character(data$dates[start_:j])

model_parameters <- list(sim           = sim,
                         burnin        = ceiling(0.5*sim),
                         N.days.fixed  =  3,
                         quantile      = c(0.1,0.9))

prior_list <- list(mu_beta    = c(0,0,0),
                   Sigma_beta = 1/2*diag(3),
                   a_sigma    = c(3,3),
                   b_sigma    = c(5/2,5/2))

result <- model(new_cases, model_parameters, prior_list)