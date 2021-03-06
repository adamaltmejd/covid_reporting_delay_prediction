##
#  building the benchmark prediction parameters for the full model
#  need to run Beta_GP_benchmark.R
##
library(foreach)
library(tidyr)
library(doParallel)
source(file.path("src", "util.r"))
source(file.path("src", "MH.R"))
source(file.path("src", "functions.R"))
source(file.path("src", "GPutil.R"))
source(file.path("src", "modelbenchmarkutil_poisson3.R"))
nclust <- 6
MCMC_sim <- 100
burnin_p = 0.5
deaths_sim <- 5
maxusage.day = 20 #must be less then N
unique.days  = 2
true.day = 8
start.predict.day = 15#16 # more then unique days

result <- readRDS(file.path("data", "processed", "processed_data.rds"))
N_T <- dim(result$detected)[1]
cl <- parallel::makeCluster(nclust,setup_strategy = "sequential", outfile="")
doParallel::registerDoParallel(cl)
foreach(j = start.predict.day:N_T)  %dopar% {
    cat('j=',j,", ", file=stdout())
    result <- readRDS(file.path("data", "processed", "processed_data.rds"))
  res_save <-benchmark_BetaGP_j(j,
                               result,true.day,
                               maxusage.day,
                               unique.days,
                               MCMC_sim,
                               burnin_p,
                               deaths_sim = deaths_sim)
  plot( colMeans(res_save$Death_est),ylim=c(20,120))
  points(apply(result$detected[1:length(colMeans(res_save$Death_est)),],1,max,na.rm=T),col='red')
  lines(apply( res_save$Death_est,2,quantile,0.95),col='green')
  lines(apply( res_save$Death_est,2,quantile,0.05),col='green')
  save(res_save,
       file = file.path("data", "simulation_results_model2", paste0("param_", result$dates_report[j], ".rds")))
}

parallel::stopCluster(cl)
