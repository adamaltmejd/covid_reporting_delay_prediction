##
#  creates preidction
#  (theta    - parameter of the detection probabilites
#  Death_est - estimated through MCMC)
#  needed for the benchmark
#  need to run Beta_GP_benchmark.R
#
#  load file:
#  data/processed/processed_data.rds (buildData.R)

#  generates the files:
#  data/simulation_results/param_lagged_YYYY-MM-DD.rdx
##

library(tidyr)
source(file.path("src", "util","util.r"))
source(file.path("src", "util","MH.R"))
source(file.path("src", "util","functions.R"))
source(file.path("src", "util","GPutil.R"))
source(file.path("src","util","MLbeta.R"))


lag <- 1
npar <- 2
MCMC_sim <- 150
#library(rGIG)
result <- readRDS(file.path("data", "processed", "processed_data.rds"))

j <- dim(result$detected)[1]
start_ = j - 20
result_j <- result

result_j$detected     <- result_j$detected[start_:j,start_:j]
result_j$dates        <- result_j$dates[start_:j]
result_j$dates_report <-  result_j$dates_report[start_:j]
result_j$dates_not_reported <- result_j$dates_not_reported[start_:j]
# remove put in deafult
alpha_Beta <-  ML_betaBin(result_j,
                   lag,
                   npar)

alpha_Beta_post <- ML_betaBin_post(result_j,
                                    lag,
                                    30,
                                    npar)
# GP prior
theta_prior <- c(0, 1)
res_save <- BetaGP_lag_fast(result_j,
                            alpha_Beta,
                            lag,
                            npar,
                            theta_prior)
