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
nclust <- 6
report.day = 14 # how long after death we wish to predicit
path.to.files <- file.path("data")
#files <- list.files('/Users/jonaswallin/Dropbox/temp/simulation_results_beta',
#                    pattern = "^param",
#                    full.names = TRUE)
files <- list.files(paste(path.to.files,"/simulation_results_old",sep=""),
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
  Date <- as.Date(result$dates_report[dim(res_save$Death_est)[2]])
  cat('file= ',format(Date) , '\n',file=stdout())

  result <- readRDS(file.path("data", "processed", "processed_data.rds"))

  pred_time = which(result$dates_report==str_extract(files[i],'[0-9]{4}-[0-9]{2}-[0-9]{2}'))
  cat('pred_time=',pred_time,"\n", file=stdout())
  N <- dim(result$detected)[1]
  report_dates <-c(result$dates_report[1:pred_time],
                   result$dates_report[pred_time] + 1:report.day)

  X_T <- setup_data(pred_time +report.day, maxusage.day, report_dates, unique.days)
  Reported <- result$detected[1:pred_time,1:pred_time]
  Reported_fill <- cbind(Reported, matrix(NA, nrow=pred_time,ncol = report.day))

  p_T <- dim(res_save$Thetas)[2]/2
  sim <- dim(res_save$Death_est)[1]



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
    Alpha_T <- matrix(NA, pred_time +report.day , pred_time +report.day)
    Beta_T  <- matrix(NA, pred_time +report.day, pred_time +report.day)

    mu <-  1/(1+exp(-X_T%*%res_save$Thetas[k,1:p_T]))
    M  <- exp(X_T%*%res_save$Thetas[k,p_T+(1:p_T)])
    Alpha_T[upper.tri(Alpha_T,diag=T)] <- M*mu
    Beta_T[upper.tri(Beta_T,diag=T)]   <- (1 - mu)*M

    Alpha_T <- Alpha_T[1:pred_time,1:(pred_time +report.day)]
    Beta_T  <- Beta_T[1:pred_time, 1:(pred_time +report.day)]


    Reported_sample <-fill.ReportBB(res_save$Death_est[k,],
                                    Alpha_T,
                                    Beta_T,
                                    Reported_fill,
                                    maxusage.day = 0)
  simulation.data$yhat[1:length(Death.date) + (length(Death.date)* (k-1))] =
    Reported_sample[row(Reported_sample)  == col(Reported_sample) -report.day
                             &  col(Reported_sample) >= pred_time]

  }

  simulation.data
}
parallel::stopCluster(cl)
states_u <- sort(unique(sim.data$state))
sim.data$yhat_mod <- sim.data$yhat
for(k in 7:length(states_u)){
  for(d in 0:14){
    data <- sim.data[     sim.data$state >= states_u[k]-10 &
                          sim.data$state < states_u[k] &
                         sim.data$days.left ==d &
                         is.na(sim.data$y) == F &
                         is.na(sim.data$yhat) == F,]
    cat('k:',k,' d:',d,' m:', median(data$y-data$yhat, na.rm=T),' mean:', mean(data$y-data$yhat, na.rm=T),'\n')
    index <- sim.data$state==states_u[k] &
            sim.data$days.left ==d
    sim.data$yhat_mod[index] <- sim.data$yhat[index] + median(data$y-data$yhat, na.rm=T)
  }

}

write_fst(sim.data, file.path("data", "processed", "model_refit.fst"))
print(sim.data[,.(median(yhat, na.rm=T),median(yhat_mod, na.rm=T)), by = .(days.left)])
SCPRS <- c()
for(k in 1:length(states_u)){
  data_k <- sim.data[sim.data$state==states_u[k],]
  dates_u <- unique(data_k$date)
 for(j in 1:length(dates_u)){
    data_kj <- data_k[data_k$date==dates_u[j],]
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
print(SCPRS[, .(mean(V2, na.rm=T),mean(V3, na.rm=T)), by = .(days.left)])
print(SCPRS[order(days.left), .(mean(EabsXX, na.rm=T), mean(V5, na.rm=T),mean(V6, na.rm=T)), by = .(days.left)])

print(SCPRS[, .(mean(V7, na.rm=T)-mean(V8, na.rm=T)), by = .(days.left)])
