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

result <- readRDS(file.path("data", "processed", "processed_data.rds"))
Reported_T = result$detected
N_T <- dim(Reported_T)[1]

deaths_est_T <- apply(Reported_T, 1, max, na.rm=T)

data_T <- newDeaths(deaths_est_T,
                   Reported_T,
                   maxusage.day =maxusage.day)
X_T <- setup_data(N_T, maxusage.day, result$dates_report, unique.days)
predicition.list <-list()


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

  save(res_save,
       file = file.path("data", "simulation_results", paste0("param_", result$dates_report[j], ".rds")))
}

parallel::stopCluster(cl)
