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
MCMC_sim <- 10000
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
                                    0,
                                    npar)
# GP prior
theta_prior <- c(3, 1)
res_post <- BetaGP_lag_post_fast(result_j,
                                 alpha_Beta_post,
                                 lag,
                                 npar,
                                 theta_prior,
                                 MCMC_sim =MCMC_sim)
theta_prior <- c(1.5, 1)
res_save <- BetaGP_lag_fast(result_j,
                            alpha_Beta,
                            lag,
                            npar,
                            theta_prior,
                            MCMC_sim  =MCMC_sim)

res_save$Death_est+res_post$Death_est
plot(result_j$dates,result_j$detected[,21])
lines(result_j$dates,apply(res_save$Death_est+res_post$Death_est,2,quantile,probs=c(0.5)))
Death_Q<-apply(res_save$Death_est+res_post$Death_est,2,quantile,probs=c(0.25,0.5,0.75))
print(rbind(result_j$detected[,21],apply(res_save$Death_est+res_post$Death_est,2,quantile,probs=c(0.25,0.5,0.75))))
