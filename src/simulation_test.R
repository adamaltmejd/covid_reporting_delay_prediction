##
#  building the benchmark prediction parameters for the full model
#  need to run Beta_GP_benchmark.R
##
set.seed(1)
library(tidyr)
library(doParallel)
source(file.path("src", "util.r"))
source(file.path("src", "MH.R"))
source(file.path("src", "functions.R"))
source(file.path("src", "GPutil.R"))
source(file.path("src", "modelbenchmarkutil_poisson.R"))
path.to.files <- file.path("data")
MCMC_sim <- 2000
burnin_p = 0.5
deaths_sim <- 5
maxusage.day = 20 #must be less then N
unique.days  = 2
true.day = 8
i = 30#16 # more then unique days
result <- readRDS(file.path("data", "processed", "processed_data.rds"))
deaths <- apply(result$detected[1:j,], 1, max, na.rm=T)
files <- list.files(paste(path.to.files,"/simulation_results_model1",sep=""),
                    pattern = "^param",
                    full.names = TRUE)
load(files[i]) #res_save_lag
j <- dim(res_save$Death_est)[2]
X_T <- setup_data(j, maxusage.day, result$dates_report[1:j], 2)
Reported <- result$detected[1:j,1:j]
data_T <- newDeaths(deaths,
                   Reported,
                   maxusage.day =maxusage.day)
p_T <- dim(res_save$Thetas)[2]/2
Alpha_T <- matrix(NA, j , j)
Beta_T  <- matrix(NA, j, j)

mu <-  1/(1+exp(-X_T%*%res_save$Thetas[k,1:p_T]))
M  <- exp(X_T%*%res_save$Thetas[k,p_T+(1:p_T)])
Alpha_T[upper.tri(Alpha_T,diag=T)] <- M*mu
Beta_T[upper.tri(Beta_T,diag=T)]   <- (1 - mu)*M

Mean_ <- matrix(0, j, j)
Mean_2 <- matrix(0, j, j)
mu_T <- Alpha_T/(Alpha_T+Beta_T)
for(i in 1:500){
report_sim <-fill.ReportBB(deaths,
                                Alpha_T,
                                Beta_T,
                       NA*Reported,
                                maxusage.day = 0)
data_ <- newDeaths(deaths,
                   report_sim,
                   maxusage.day =maxusage.day)
Mean_ <- Mean_ + data_$report.new
Mean_2 <- Mean_2 + (data_$report.new)^2
}
Mean_ <- Mean_/500
Mean_2 <- Mean_2/500
V_2    <- sqrt(Mean_2 - Mean_^2)
sst <- (Mean_-data_T$report.new)/V_2
result_sim <- result
result_sim$detected <- report_sim
result_sim$dates <- result$dates[1:j]
result_sim$dates_report <- result$dates_report[1:j]
res_sim <-benchmark_BetaGP_j(j,
                              result_sim,true.day,
                              maxusage.day,
                              unique.days,
                              MCMC_sim,
                              burnin_p,
                              deaths_sim = deaths_sim)
plot(colMeans(res_sim$Death_est), ylim=c(10,120))
points(deaths, col='red')
points(colMeans(res_save$Death_est), col='blue', pch=1)
