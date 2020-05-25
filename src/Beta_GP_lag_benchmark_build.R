##
#  building the benchmark prediction parameters for the full model
#  need to run Beta_GP_benchmark.R
##
library(foreach)
library(doParallel)
source(file.path("src", "util.r"))
source(file.path("src", "MH.R"))
source(file.path("src", "functions.R"))
source(file.path("src", "GPutil.R"))
source(file.path("src", "MLbeta.R"))
nclust <- 6
lag <- 3
MCMC_sim <- 15000
burnin_p = 0.5
deaths_sim <- 5
start.predict.day = 16#16 # more then unique days
predicition.list <-list()


cl <- parallel::makeCluster(nclust,setup_strategy = "sequential", outfile="")
doParallel::registerDoParallel(cl)
foreach(j = start.predict.day:N_T)  %dopar% {
    cat('j=',j,", ", file=stdout())
    result <- readRDS(file.path("data", "processed", "processed_data.rds"))
  res_save<-benchmark_BetaGP_lag_j(j,
                                    result = result,
                                    lag   = lag,
                                    MCMC_sim = MCMC_sim,
                                    burnin_p = burnin_p,
                                    deaths_sim = deaths_sim,
                                    prior = c(0))

  save(res_save,
       file = file.path("data", "simulation_results", paste0("param_lagged_", result$dates_report[j], ".rds")))
}

parallel::stopCluster(cl)
