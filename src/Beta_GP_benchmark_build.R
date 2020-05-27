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
source(file.path("src", "modelbenchmarkutil_poisson.R"))
nclust <- 6
MCMC_sim <- 15000
burnin_p = 0.5
deaths_sim <- 5
maxusage.day = 20 #must be less then N
unique.days  = 5
true.day = 5
start.predict.day = 40#16 # more then unique days



cl <- parallel::makeCluster(nclust,setup_strategy = "sequential", outfile="")
doParallel::registerDoParallel(cl)
foreach(j = start.predict.day:N_T)  %dopar% {
    cat('j=',j,", ", file=stdout())
    result <- readRDS(file.path("data", "processed", "processed_data.rds"))
  res_save<-benchmark_BetaGP_j(j,
                               result,true.day,
                               maxusage.day,
                               unique.days,
                               MCMC_sim,
                               burnin_p,
                               deaths_sim = deaths_sim)

  plot(colMeans(exp(res_save$theta_GP)),ylim=c(0,120))
  points(apply(replace_na(result$detected[1:j,1:j],0),1,max, na.rm=T),col='green')
  points(colMeans(res_save$Death_est),col='blue')
  lines(apply(res_save$Death_est,2,quantile, c(0.05,0.95))[1,])
  lines(apply(res_save$Death_est,2,quantile, c(0.05,0.95))[2,])
  points(apply(result$detected[1:j,],1,max, na.rm=T),col='red')
  save(res_save,
       file = file.path("data", "simulation_results", paste0("param_", result$dates_report[j], ".rds")))
}

parallel::stopCluster(cl)
