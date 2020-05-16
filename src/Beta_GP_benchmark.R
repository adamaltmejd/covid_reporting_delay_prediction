###
# benchmarking 
# 
###
source('util.r')
path.to.files <- file.path("..","data")
files <- list.files(paste(path.to.files,"/simulation_result",sep=""),
                             pattern = "^param",
                             full.names = TRUE)
for(i in 1:length(files)){

  load(files[i])
  load(paste(path.to.files,"/result.RData",sep=""))
  
  Reported_T = result$detected 
  pred_time = dim(Death_est)[2]
  N_T <- dim(Reported_T)[1]
  if(pred_time < N_T){
    deaths_est_T <- apply(Reported_T, 1, max, na.rm=T)
    data_T <- newDeaths(deaths_est_T,
                        Reported_T, 
                        maxusage.day =maxusage.day)
    X_T <- setup_data(N_T, maxusage.day, result$dates_report, unique.days)
    Alpha_T <- matrix(NA, N_T,N_T)
    Beta_T  <- matrix(NA, N_T,N_T)
    
    
    
    
    Reported <- result$detected[1:pred_time,1:pred_time]
    Reported_fill <- cbind(Reported, matrix(NA, nrow=pred_time,ncol = (N_T-pred_time)))
    
    p <- dim(Thetas)[2]/2
    sim <- dim(Thetas)[2]
    for(k in 1:sim)
      beta_1 <- Thetas[k,1:p]
      beta_2 <- Thetas[k,(p+1):(2*p)]
      Alpha_T[upper.tri(Alpha_T,diag=T)] <- exp(X_T%*%beta_1)
      Beta_T[upper.tri(Beta_T,diag=T)]   <- exp(X_T%*%beta_2)
      Alpha_T <- Alpha_T[1:pred_time,1:N_T]
      Beta_T  <- Beta_T[1:pred_time, 1:N_T]
      Reported_sample <-fill.ReportBB(Death_est[k,], 
                                      Alpha_T,
                                      Beta_T, 
                                      Reported_fill, 
                                      maxusage.day = true.day)
      #pred_set[,,i-burnin + 1] <- Reported_sample[,(j+1):dim(Reported_sample)[2]]
    }
}
  
 