###
# benchmarking 
# 
###
library(foreach)
source('util.r')
nclust <- 4
path.to.files <- file.path("..","data")
files <- list.files(paste(path.to.files,"/simulation_result",sep=""),
                             pattern = "^param",
                             full.names = TRUE)
alpha_CI <- 0.1
load(paste(path.to.files,"/result.RData",sep=""))
cl <- parallel::makeCluster(nclust)
doParallel::registerDoParallel(cl)
index_files <-1:(length(files)-1)
#for(i in index_files){
result_par <-foreach(i = index_files)  %dopar% {
  library(Matrix)

  load(files[i])
  load(paste(path.to.files,"/result.RData",sep=""))
  
  Reported_T = result$detected 
  pred_time = dim(Death_est)[2]
  N_T <- dim(Reported_T)[1]
  deaths_est_T <- apply(Reported_T, 1, max, na.rm=T)
  data_T <- newDeaths(deaths_est_T,
                      Reported_T, 
                      maxusage.day =maxusage.day)
  X_T <- setup_data(N_T, maxusage.day, result$dates_report, unique.days)
  
  
  
  
  Reported <- result$detected[1:pred_time,1:pred_time]
  Reported_fill <- cbind(Reported, matrix(NA, nrow=pred_time,ncol = (N_T-pred_time)))
  
  p <- dim(Thetas)[2]/2
  sim <- dim(Thetas)[1]
  pred_set <-array(NA, dim = c(pred_time,N_T-pred_time,sim)) 
  for(k in 1:sim){
    beta_1 <- Thetas[k,1:p]
    beta_2 <- Thetas[k,(p+1):(2*p)]
    Alpha_T <- matrix(NA, N_T,N_T)
    Beta_T  <- matrix(NA, N_T,N_T)
    
    Alpha_T[upper.tri(Alpha_T,diag=T)] <- exp(X_T%*%beta_1)
    Beta_T[upper.tri(Beta_T,diag=T)]   <- exp(X_T%*%beta_2)
    Alpha_T <- Alpha_T[1:pred_time,1:N_T]
    Beta_T  <- Beta_T[1:pred_time, 1:N_T]
    Reported_sample <-fill.ReportBB(Death_est[k,], 
                                    Alpha_T,
                                    Beta_T, 
                                    Reported_fill, 
                                    maxusage.day = true.day)
    pred_set[,,k] <- Reported_sample[,(pred_time+1):N_T]
  }
  
  CI <- apply(pred_set,c(1,2), quantile,probs=c(alpha_CI/2,0.5,1-alpha_CI/2))
  CI_low <- as.matrix(CI[1,,])
  colnames( CI_low)<- format(result$dates_report[(pred_time+1):N_T],"%Y-%m-%d")
  rownames( CI_low)<- format(result$dates[1:pred_time],"%Y-%m-%d")
  
  med <- as.matrix(CI[2,,])
  colnames( med)<- format(result$dates_report[(pred_time+1):N_T],"%Y-%m-%d")
  rownames( med)<- format(result$dates[1:pred_time],"%Y-%m-%d")
  
  CI_up <- as.matrix(CI[3,,])
  colnames( CI_up)<- format(result$dates_report[(pred_time+1):N_T],"%Y-%m-%d")
  rownames( CI_up)<- format(result$dates[1:pred_time],"%Y-%m-%d")
  index <- format(result$dates_report[pred_time],"%Y-%m-%d")
  predicition_by_reported_day <-list()
  predicition_by_reported_day$CI_low <- CI_low
  predicition_by_reported_day$CI_up <- CI_up
  predicition_by_reported_day$med <- med
  predicition_by_reported_day$Truth <- as.matrix(result$detected[1:pred_time,(pred_time+1):N_T])
  
  colnames( predicition_by_reported_day$Truth)<- format(result$dates_report[(pred_time+1):N_T],"%Y-%m-%d")
  rownames( predicition_by_reported_day$Truth)<- format(result$dates[1:pred_time],"%Y-%m-%d")
  SCPRS <- matrix(nrow=pred_time, ncol = N_T-pred_time)
  for(k_death in 1:pred_time){
    for(k_report in 1:(N_T-pred_time)){
      
      predSample <-pred_set[k_death,k_report,]
      P1 <- predSample[1:ceiling(sim/3)]
      P2 <- predSample[(ceiling(sim/3) + 1):(floor(2*sim/3) )]
      P3 <- predSample[(floor(2*sim/3) + 1):(floor(2*sim/3) + length(P2))]
      EabsXX <- mean(abs(sample(P2)-sample(P3)))
      Y <- predicition_by_reported_day$Truth[k_death,k_report]
      EabsYX <- mean(abs(Y-P1))
      SCPRS[k_death, k_report] <- -EabsYX/EabsXX - 0.5*log(EabsXX) 
    }
  }
  colnames( SCPRS)<- format(result$dates_report[(pred_time+1):N_T],"%Y-%m-%d")
  rownames( SCPRS)<- format(result$dates[1:pred_time],"%Y-%m-%d")
  predicition_by_reported_day$SCPRS <- SCPRS
  predicition_by_reported_day
}

parallel::stopCluster(cl)
names(result_par) <- result$dates_report[start.predict.day+index_files-1]
save(
     result_par,
     file=paste(path.to.files,"/simulation_result/prediction_valdiation.RData",sep="")
     )
date_predict <- '2020-05-03'
date_death   <- '2020-05-01'
rbind(result_par[[date_predict]]$CI_low[date_death,],
      result_par[[date_predict]]$med[date_death,],
      result_par[[date_predict]]$Truth[date_death,],
      result_par[[date_predict]]$CI_up[date_death,],
      result_par[[date_predict]]$SCPRS[date_death,])

 