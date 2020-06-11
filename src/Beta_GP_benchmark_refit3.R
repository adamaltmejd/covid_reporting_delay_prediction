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


  X_mixed <- setup_data_mixed_effect(pred_time +report.day,
                                     2,
                                     result$dates_report[1:(pred_time +report.day)])
  X_mu <- setup_data(pred_time +report.day,
                     maxusage.day,
                     result$dates_report[1:(pred_time +report.day)], 3)
  X_M  <- X_mu[,1:2]
  X_longlag <- rowSums(X_mixed)>0
  X_trend <- setup_all_days(pred_time +report.day)/7
  #X_mu <- cbind(X_mu,)
  X_M  <- cbind(X_M[,1],X_M[,2], X_longlag)

  p <- dim(X_mu)[2]
  p1 <- dim(X_M)[2]
  p2 <- dim(X_mixed)[2]


  Reported <- result$detected[1:pred_time,1:pred_time]
  Reported_fill <- cbind(Reported, matrix(NA, nrow=pred_time,ncol = report.day))

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
    beta_mu    <- res_save$Thetas[k,1:p]
    beta_M     <- res_save$Thetas[k,p + (1:p1)]
    mu <- 1/(1+exp(-X_mu%*%beta_mu ))
    M  <- exp(X_M%*%beta_M)
    Alpha_T <- matrix(NA, pred_time +report.day , pred_time +report.day)
    Beta_T  <- matrix(NA, pred_time +report.day, pred_time +report.day)

    Alpha_T[upper.tri(Alpha_T,diag=T)] <- M*mu
    Beta_T[upper.tri(Beta_T,diag=T)]   <- (1 - mu)*M

    Alpha_T<- Alpha_T[1:pred_time,1:(pred_time +report.day)]
    Beta_T  <- Beta_T[1:pred_time, 1:(pred_time +report.day)]


    Reported_sample <-fill.ReportBB(res_save$Death_est[k,],
                                    Alpha_T,
                                    Beta_T,
                                    Reported_fill,
                                    maxusage.day = 0)
  simulation.data$yhat[1:length(Death.date) + (length(Death.date)* (k-1))] =
    Reported_sample[row(Reported_sample)  == col(Reported_sample) -report.day
                             &  col(Reported_sample) >= pred_time]
  simulation.data$Nhat[1:length(Death.date) + (length(Death.date)* (k-1))] =res_save$Death_est[k,]

  }

  simulation.data
}
parallel::stopCluster(cl)
states_u <- sort(unique(sim.data$state))
sim.data$yhat_mod <- sim.data$yhat
sim.data$yhat_mod2 <- sim.data$yhat
for(k in 7:length(states_u)){
  for(d in 0:14){
    data <- sim.data[     sim.data$state >= states_u[k]-10 -(14-d) &
                          sim.data$state < states_u[k] - (14-d) &
                         sim.data$days.left ==d &
                         is.na(sim.data$y) == F &
                         is.na(sim.data$yhat) == F,]
    cat('k:',k,' d:',d,' m:', median(data$y-data$yhat, na.rm=T),' mean:', mean(data$y-data$yhat, na.rm=T),'\n')
    index <- sim.data$state==states_u[k] &
            sim.data$days.left ==d
    if(is.na( median(data$y-data$yhat, na.rm=T))==F & d < 10){
      sim.data$yhat_mod[index] <- sim.data$yhat[index] + median(data$y-data$yhat, na.rm=T)
      dats <- unique(data$date)
      data$yhat2 <- NA
      for(dat in 1:length(dats)){
        index <- data$date==dats[dat]
        data$yhat2[index] <- sample(data$yhat[index])
      }
      sim.data$yhat_mod2 <- sim.data$yhat_mod *sd(data$y-data$yhat, na.rm=T)/sd(data$yhat2-data$yhat, na.rm=T)
    }

  }

}

write_fst(sim.data, file.path("data", "processed", "model1_refit.fst"))
model_benchmark       <- c()
model_benchmark_bias  <- c()
SCPRS                 <- c()
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
    l <- min(length(P1),length(P3))
    EabsXX <- mean(abs(sample(P1[1:l])-sample(P3[1:l])))
    if(EabsXX == 0){
      EabsXX <- 1
    }
    EabsYX <- mean(abs(data_kj$y[1]-P2))
    CI_1 <- as.vector(quantile(data_kj$yhat,prob=c(0.05,0.95),na.rm=T))
    M1   <- median(data_kj$yhat, na.rm=T)


    P1_mod <- data_kj$yhat_mod[1:ceiling(sim/3)]
    P2_mod <- data_kj$yhat_mod[(ceiling(sim/3) + 1):(floor(2*sim/3) )]
    P3_mod <- data_kj$yhat_mod[(floor(2*sim/3) + 1):(floor(2*sim/3) + length(P2))]
    EabsXX_mod <- mean(abs(sample(P1_mod[1:l])-sample(P3_mod[1:l])))

    if(is.na(EabsXX_mod)){
      EabsYX_mod <- EabsYX
      EabsXX_mod <- EabsXX
      CI_2 <- CI_1
      M2 <- M1
    }else{
      if(EabsXX_mod == 0){
        EabsXX_mod <- 1
      }
      EabsYX_mod <- mean(abs(data_kj$y[1]-P2_mod))

      CI_2 <- as.vector(quantile(data_kj$yhat_mod,prob=c(0.05,0.95),na.rm=T))
      M2   <- median(data_kj$yhat_mod, na.rm=T)
    }

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
model_benchmark_bias$days_left =  model_benchmark_bias$date -model_benchmark_bias$state +14

names(model_benchmark) <- c('state','date','days_left','target','SCRPS','predicted_deaths','ci_lower','ci_upper')
model_benchmark<-model_benchmark[,c('state','date','target','days_left','ci_upper','ci_lower','predicted_deaths','SCRPS')]
model_benchmark$days_left =  model_benchmark$date -model_benchmark$state +14

write_fst(model_benchmark, file.path("data", "processed", "model_benchmark.fst"))
write_fst(model_benchmark_bias, file.path("data", "processed", "model_benchmark.fst"))

sim_today <- sim.data[state==max(states_u),]
alpha=0.1
sim_today<- sim_today[,.(quantile(yhat_mod,prob=c(1-alpha/2), na.rm=T),
             quantile(yhat_mod,prob=c(alpha/2), na.rm=T),
             quantile(yhat_mod,prob=c(0.5), na.rm=T)), by = .(date)]
names(sim_today) <- c("date","ci_upper","ci_lower","predicted_deaths")
write_fst(sim_today, file.path("data", "processed", "model_latest.fst"))


