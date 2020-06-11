##
#  building the benchmark prediction parameters for the full model
#  need to run Beta_GP_benchmark.R
##
library(foreach)
library(doParallel)
library(tidyr)
source(file.path("src", "util.r"))
source(file.path("src", "MH.R"))
source(file.path("src", "functions.R"))
source(file.path("src", "GPutil.R"))
source(file.path("src","lagged" ,"MLbeta.R"))
nclust <- 6
lag <- 1
MCMC_sim <- 15000
burnin_p = 0.3
deaths_sim <- 5
true.day <- 5
maxusage.day <- 20
start.predict.day = 20#16 # more then unique days

predicition.list <-list()
result <- readRDS(file.path("data", "processed", "processed_data.rds"))

cl <- parallel::makeCluster(nclust,setup_strategy = "sequential", outfile="")
doParallel::registerDoParallel(cl)
foreach(j = start.predict.day:dim(result$detected)[1], .errorhandling='remove')  %dopar% {
    library(rGIG)
    cat('j=',j,", ", file=stdout())
    result <- readRDS(file.path("data", "processed", "processed_data.rds"))
    res_save_post_lag<-benchmark_BetaGP_post_lag_j(j,
                                                   result,
                                                   lag,
                                                   true.day,
                                                   maxusage.day,
                                                   MCMC_sim,
                                                   burnin_p,
                                                   deaths_sim,
                                                   prior = 0,
                                                   npars = 2)
  res_save_lag<-benchmark_BetaGP_lag_j(j,
                                      result = result,
                                      lag   = lag,
                                      MCMC_sim = MCMC_sim,
                                      burnin_p = burnin_p,
                                      deaths_sim = deaths_sim,
                                      prior = 0,
                                      npar = 2)
  save(res_save_lag,res_save_post_lag,
       file = file.path("data", "simulation_results", paste0("param_lagged_", result$dates_report[j], ".rds")))

  if(0){
    j <- res_save_lag$j
  result_j <- result
  result_j$detected     <- result_j$detected[1:j,1:j]
  result_j$dates        <- result_j$dates[1:j]
  result_j$dates_report <-  result_j$dates_report[1:j]
  report_j <- splitlag(result_j$detected,as.Date(result_j$dates_report), lag)
  report   <- splitlag(result$detected,as.Date(result$dates_report), lag)
  plot(colMeans(exp(t(res_save_lag$A%*%t(res_save_lag$theta_GP)))),ylim=c(0,100))
  points(apply(report_j$Reported_T[1:j,],1,max, na.rm=T),col='green')
  points(colMeans(res_save_lag$Death_est),col='blue')
  lines(apply(res_save_lag$Death_est,2,quantile, c(0.05,0.95))[1,])
  lines(apply(res_save_lag$Death_est,2,quantile, c(0.05,0.95))[2,])
  points(apply(report$Reported_T[1:j,],1,max, na.rm=T),col='red')


  plot(colMeans(exp(t(res_save_post_lag$A%*%t(res_save_post_lag$theta_GP)))),ylim=c(0,100))
  points(apply(report_j$Reported_O[1:j,],1,max, na.rm=T),col='green')
  points(colMeans(res_save_post_lag$Death_est),col='blue')
  lines(apply(res_save_post_lag$Death_est,2,quantile, c(0.05,0.95))[1,])
  lines(apply(res_save_post_lag$Death_est,2,quantile, c(0.05,0.95))[2,])
  points(apply(report$Reported_O[1:j,],1,max, na.rm=T),col='red')

  plot(colMeans(exp(t(res_save_post_lag$A%*%t(res_save_post_lag$theta_GP))))+
         colMeans(exp(t(res_save_lag$A%*%t(res_save_lag$theta_GP))))
       ,ylim=c(0,120))
  points(apply(replace_na(report_j$Reported_O[1:j,1:j],0),1,max)+
           apply(report_j$Reported_T[1:j,],1,max, na.rm=T),col='green')
  points(colMeans(res_save_post_lag$Death_est)+
           colMeans(res_save_lag$Death_est),col='blue')
  lines(apply(res_save_post_lag$Death_est + res_save_lag$Death_est,2,quantile, c(0.05,0.95))[1,])
  lines(apply(res_save_post_lag$Death_est + res_save_lag$Death_est,2,quantile, c(0.05,0.95))[2,])
  points(apply(result$detected[1:j,],1,max, na.rm=T),col='red')
  }
}

parallel::stopCluster(cl)
