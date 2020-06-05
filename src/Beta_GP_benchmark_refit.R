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
maxusage.day = 20
unique.days  = 2
true.day = 5
path.to.files <- file.path("data")
#files <- list.files('/Users/jonaswallin/Dropbox/temp/simulation_results_beta',
#                    pattern = "^param",
#                    full.names = TRUE)
files <- list.files(paste(path.to.files,"/simulation_results_model1",sep=""),
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
                                yhat       = NA,
                                Nhat       = NA)



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
    simulation.data$Nhat[1:length(Death.date) + (length(Death.date)* (k-1))] =res_save$Death_est[k,]
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
    data <- sim.data[    sim.data$state >= states_u[k]-10 -(14-d) &
                         sim.data$state < states_u[k] - (14-d) &
                         sim.data$days.left ==d &
                         is.na(sim.data$y) == F &
                         is.na(sim.data$yhat) == F,]
    cat('k:',k,' d:',d,' m:', median(data$y-data$yhat, na.rm=T),' mean:', mean(data$y-data$yhat, na.rm=T),'\n')
    index <- sim.data$state==states_u[k] &
            sim.data$days.left ==d
    if(is.na(median(data$y-data$yhat, na.rm=T))==F){
      sim.data$yhat_mod[index] <- sim.data$yhat[index] + median(data$y-data$yhat, na.rm=T)
    }

  }

}

write_fst(sim.data, file.path("data", "processed", "model_refit.fst"))
print(sim.data[,.(median(yhat, na.rm=T),median(yhat_mod, na.rm=T)), by = .(days.left)])
SCPRS <- c()
model_benchmark_bias <- c()
model_benchmark      <- c()
for(k in 1:length(states_u)){
  data_k <- sim.data[sim.data$state==states_u[k],]
  dates_u <- unique(data_k$date)
 for(j in 1:length(dates_u)){
    data_kj <- data_k[data_k$date==dates_u[j],]
    if(is.na(data_kj$y[1])==T)
      next
    sim <- length(data_kj$yhat)
    P1 <- data_kj$yhat[1:(floor(sim/3) )]
    P2 <- data_kj$yhat[(floor(sim/3) + 1):(floor(2*sim/3) )]
    P3 <- data_kj$yhat[(floor(2*sim/3) + 1):(floor(2*sim/3) + length(P1))]
    EabsXX <- mean(abs(sample(P1)-sample(P3)))
    if(EabsXX == 0){
      EabsXX <- 1
    }
    EabsYX <- mean(abs(data_kj$y[1]-P2))

    P1_mod <- data_kj$yhat_mod[1:(floor(sim/3) )]
    P2_mod <- data_kj$yhat_mod[(floor(sim/3) + 1):(floor(2*sim/3) )]
    P3_mod <- data_kj$yhat_mod[(floor(2*sim/3) + 1):(floor(2*sim/3) + length(P1_mod))]

    EabsXX_mod <- mean(abs(sample(P1_mod)-sample(P3_mod)))
    I <- -1-2*prod(sign(quantile(data_kj$yhat_mod,prob=c(0.05,0.95), na.rm=T)))
    if(is.na(EabsXX_mod)==T){
      EabsXX_mod <- EabsXX
      EabsYX_mod <- EabsYX
      I_mod <- I
    }else{
      if(EabsXX_mod == 0){
        EabsXX_mod <- 1
      }
      EabsYX_mod <- mean(abs(data_kj$y[1]-P2_mod))
      I_mod<- prod(sign(quantile(data_kj$yhat_mod,prob=c(0.05,0.95), na.rm=T)))
    }
    CI_1 <- as.vector(quantile(data_kj$yhat,prob=c(0.05,0.95),na.rm=T))
    CI_2 <- as.vector(quantile(data_kj$yhat_mod,prob=c(0.05,0.95),na.rm=T))
    M1   <- median(data_kj$yhat, na.rm=T)
    M2   <- median(data_kj$yhat_mod, na.rm=T)
    SCPRS <- rbind(SCPRS,
                   cbind(data_kj[1,1:4],
                         -EabsYX/EabsXX - 0.5*log(EabsXX),
                         -EabsYX_mod/EabsXX_mod - 0.5*log(EabsXX_mod),
                         CI_1[1],CI_1[2],CI_2[1],CI_2[2],
                         (CI_1[1]-data_kj$y[1])*(CI_1[2]-data_kj$y[1])<=0,
                          (CI_2[1]-data_kj$y[1])*(CI_2[2]-data_kj$y[1])<=0,
                   M1,M2))

    model_benchmark      <-rbind(model_benchmark,
                                 cbind(data_kj[1,1:4],-EabsYX/EabsXX - 0.5*log(EabsXX),
                                 M1,
                                 CI_1[1],
                                 CI_1[2]))
    model_benchmark_bias <-rbind(model_benchmark_bias,
                                 cbind(data_kj[1,1:4],
                                       -EabsYX_mod/EabsXX_mod - 0.5*log(EabsXX_mod),
                                 M2,
                                 CI_2[1],
                                 CI_2[2]))
    }
}
names(model_benchmark_bias) <- c('state','date','days_left','target','SCRPS','predicted_deaths','ci_lower','ci_upper')
model_benchmark_bias <- model_benchmark_bias[,c('state','date','target','days_left','ci_upper','ci_lower','predicted_deaths','SCRPS')]
names(model_benchmark) <- c('state','date','days_left','target','SCRPS','predicted_deaths','ci_lower','ci_upper')
model_benchmark<-model_benchmark[,c('state','date','target','days_left','ci_upper','ci_lower','predicted_deaths','SCRPS')]

write_fst(model_benchmark, file.path("data", "processed", "model_benchmark.fst"))
write_fst(model_benchmark_bias, file.path("data", "processed", "model_benchmark.fst"))
