###
# benchmarking
#
###
library( stringr)
library(foreach)
library(tidyr)
library(fst)
library(data.table)
source("src/util.r")
nclust <- 4
report.day = 14 # how long after death we wish to predicit
start.predict.day=25
path.to.files <- file.path("data")
#files <- list.files('/Users/jonaswallin/Dropbox/temp/simulation_results_beta',
#                    pattern = "^param",
#                    full.names = TRUE)
files <- list.files(paste(path.to.files,"/simulation_results",sep=""),
                             pattern = "^param",
                             full.names = TRUE)

index_files <- 1:(length(files)-1)
result <- readRDS(file.path("data", "processed", "processed_data.rds"))


cl <- parallel::makeCluster(nclust,setup_strategy = "sequential", outfile="")
doParallel::registerDoParallel(cl)

#for(i in index_files){
sim.data <- foreach(i = index_files, .combine='rbind')  %dopar% {
  library(Matrix)
  library(tidyr)
  library(fst)
  library(data.table)
  library(stringr)
  load(files[i]) #res_save_lag

  cat('i=',i,", ", file=stdout())
  Date <- as.Date(result$dates_report[dim(res_save_lag$Death_est)[2]])
  cat('file= ',format(Date) , '\n',file=stdout())
  write.data <- file.exists(file.name)==F

  result <- readRDS(file.path("data", "processed", "processed_data.rds"))
  result$detected[14, 21:dim(result$detected)[2]] <- result$detected[14, 21:dim(result$detected)[2]] - 15

  pred_time = which(result$dates_report==str_extract(files[i],'[0-9]{4}-[0-9]{2}-[0-9]{2}'))
  cat('pred_time=',pred_time,"\n", file=stdout())
  N <- dim(result$detected)[1]
  lag <- res_save_lag$lag
  report <- splitlag(result$detected,as.Date(result$dates_report) ,res_save_lag$lag)
  report_dates <-c(result$dates_report[1:pred_time],
                   result$dates_report[pred_time] + 1:report.day)
  X_T <- setup_data_lag(pred_time +report.day, res_save_lag$npar, report_dates, res_save_lag$npar)
  X_O <- setup_data_postlag2(pred_time +report.day, res_save_lag$lag, res_save_post_lag$npars, report_dates)
  Reported <- result$detected[1:pred_time,1:pred_time]
  Reported_fill <- cbind(Reported, matrix(NA, nrow=pred_time,ncol = report.day))

  p_T <- length(res_save_lag$theta)/2
  p_O <- dim(X_O)[2]
  sim <- dim(res_save_lag$Death_est)[1]

  Alpha_T <- matrix(NA, pred_time +report.day , pred_time +report.day)
  Beta_T  <- matrix(NA, pred_time +report.day, pred_time +report.day)

  mu <-  1/(1+exp(-X_T%*%res_save_lag$theta[1:p_T]))
  M  <- exp(X_T%*%res_save_lag$theta[p_T+(1:p_T)])
  Alpha_T[upper.tri(Alpha_T,diag=T)] <- M*mu
  Beta_T[upper.tri(Beta_T,diag=T)]   <- M*(1 - mu)

  Alpha_T <- Alpha_T[1:pred_time,1:(pred_time +report.day)]
  Beta_T  <- Beta_T[1:pred_time, 1:(pred_time +report.day)]
  Alpha_O <- matrix(NA, pred_time +report.day,pred_time +report.day)
  Beta_O  <- matrix(NA, pred_time +report.day,pred_time +report.day)
  mu <-  1/(1+exp(-X_O%*%res_save_post_lag$theta[1:p_O]))
  M  <- exp(X_O%*%res_save_post_lag$theta[p_O+(1:p_O)])
  Alpha_O[upper.tri(Alpha_O,diag=T)]  <- M*mu
  Beta_O[upper.tri(Alpha_O,diag=T)]   <- M*(1- mu)
  Alpha_O <- Alpha_O[1:pred_time,1:(pred_time +report.day)]
  Beta_O   <- Beta_O[1:pred_time, 1:(pred_time +report.day)]


  state      <- as.Date(result$dates_report[pred_time])
  Death.date <- result$dates_report[max(1,pred_time-report.day):pred_time]
  y          <- result$detected[row(result$detected)  == col(result$detected) - report.day
                                &  col(result$detected) >= pred_time
                                & col(result$detected) <= pred_time + report.day]
  y <- c(y, rep(NA, max(0,length(Death.date)-length(y))))
  simulation.data <- data.table(state      = as.Date(result$dates_report[pred_time]),
                                date       = rep(Death.date,sim),
                                days.left  = rep(as.Date(result$dates_report[pred_time]) - Death.date  ,sim),
                                y          = rep(y,sim),
                                yhat       = NA)



  for(k in 1:sim){


    Reported_sample_lag <-fill.ReportBB.lag(res_save_lag$Death_est[k,],
                                    Alpha_T,
                                    Beta_T,
                                    Reported_fill,
                                    report_dates,
                                    lag = lag)
    Reported_sample_post_lag <-fill.ReportBB.postlag(res_save_post_lag$Death_est[k,],
                                            Alpha_O,
                                            Beta_O,
                                            Reported_sample_lag,
                                            report_dates,
                                            lag = lag)
  simulation.data$yhat[1:length(Death.date) + (length(Death.date)* (k-1))] =
    Reported_sample_post_lag[row(Reported_sample_post_lag)  == col(Reported_sample_post_lag) -report.day
                             &  col(Reported_sample_post_lag) >= pred_time]

  }

  simulation.data
}
parallel::stopCluster(cl)
states_u <- sort(unique(sim.data$state))
sim.data$yhat_mod <- sim.data$yhat
for(k in 5:length(states_u)){
  for(d in 0:14){
    data <- sim.data[     sim.data$state >= states_u[k]-14 &
                          sim.data$state < states_u[k] &
                         sim.data$days.left ==d &
                         is.na(sim.data$y) == F &
                         is.na(sim.data$yhat) == F,]
    cat('k:',k,' d:',d,' m:', median(data$y-data$yhat),'\n')
    index <- sim.data$state==states_u[k] &
            sim.data$days.left ==d
    sim.data$yhat_mod[index] <- sim.data$yhat[index] + mean(data$y-data$yhat)
  }

}
SCPRS <- c()
for(k in 1:length(states_u)){
  data_k <- sim.data[state==states_u[k],]
  dates_u <- unique(data_k$date)
 for(j in 1:length(dates_u)){
    data_kj <- data_k[date==dates_u[j],]
    if(is.na(data_kj$y[1])==T)
      next
    sim <- length(data_kj$yhat)
    P1 <- data_kj$yhat[1:ceiling(sim/3)]
    P2 <- data_kj$yhat[(ceiling(sim/3) + 1):(floor(2*sim/3) )]
    P3 <- data_kj$yhat[(floor(2*sim/3) + 1):(floor(2*sim/3) + length(P2))]
    EabsXX <- mean(abs(sample(P1)-sample(P3)))
    if(EabsXX == 0){
      EabsXX <- 1
    }
    EabsYX <- mean(abs(data_kj$y[1]-P2))

    P1_mod <- data_kj$yhat_mod[1:ceiling(sim/3)]
    P2_mod <- data_kj$yhat_mod[(ceiling(sim/3) + 1):(floor(2*sim/3) )]
    P3_mod <- data_kj$yhat_mod[(floor(2*sim/3) + 1):(floor(2*sim/3) + length(P2))]
    EabsXX_mod <- mean(abs(sample(P1_mod)-sample(P3_mod)))
    if(EabsXX_mod == 0){
      EabsXX_mod <- 1
    }
    EabsYX_mod <- mean(abs(data_kj$y[1]-P2_mod))

    SCPRS <- rbind(SCPRS,
                   cbind(data_kj[1,1:3],
                         -EabsYX/EabsXX - 0.5*log(EabsXX),
                         -EabsYX_mod/EabsXX_mod - 0.5*log(EabsXX_mod),
                         EabsXX,
                         -EabsYX ,
                         -EabsYX_mod,
                         median(P1_mod),
                         median(P1)))

    }
}
print(cbind(SCPRS[order(days.left,decreasing = T), .(mean(V2, na.rm=T),mean(V3, na.rm=T)), by = .(days.left)],benchmark[state>='2020-04-11', .(mean(SCRPS, na.rm=T)), by = .(days_left)]))
