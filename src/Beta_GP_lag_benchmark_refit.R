###
# benchmarking
# using parameters(theta     - parameters of the detection probabilites
#                  Death_est - true number of deaths)
# from Beta_GP_lag_benchmark_build.R to generate
# preidcition of the reporting
#
# load files:
#  data/processed/processed_data.rds   (buildData.R)
#  data/simulation_results/param_lagged_YYYY-MM-DD.rdx (Beta_GP_lag_benchmark_build.R)
#
# generates files:
# data/processed/model_latest.fst
# data/processed/model_lag_refit.fst
# data/processed/model_benchmark.fst
# data/processed/model_benchmark_bias.fst
###
library( stringr)
library(foreach)
library(tidyr)
library(fst)
library(data.table)
source(file.path("src","util","util.r"))
nclust <- 6
report.day = 14 # how long after death we wish to predicit
start.predict.day=20
path.to.files <- file.path("data")
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

  result <- readRDS(file.path("data", "processed", "processed_data.rds"))

  pred_time = which(result$dates_report==str_extract(files[i],'[0-9]{4}-[0-9]{2}-[0-9]{2}'))
  cat('pred_time=',pred_time,"\n", file=stdout())
  N <- dim(result$detected)[1]
  lag <- res_save_lag$lag
  report <- splitlag(result$detected,as.Date(result$dates_report) ,res_save_lag$lag)
  report_dates <-c(result$dates_report[1:pred_time],
                   result$dates_report[pred_time] + 1:report.day)
  X_T <- setup_data_lag(pred_time +report.day, res_save_lag$npar, report_dates, res_save_lag$npar)

  X_O <- setup_data_postlag2(pred_time +report.day,
                             res_save_lag$lag,
                             res_save_post_lag$npars,
                             report_dates)
  deaths_est <- apply(report_j$Reported_O, 1, function(x) { if(all(is.na(x))){return(0)};
    max(x,na.rm=T)})

  p_O <- dim(X_O)[2]
  X_OM  <- cbind(rep(1,dim(X_O)[1]),rowSums(X_O[,4:5])>0)
  p_OM <- dim(X_OM)[2]
  Reported <- result$detected[1:pred_time,1:pred_time]
  Reported_fill <- cbind(Reported, matrix(NA, nrow=pred_time,ncol = report.day))

  p_T <- length(res_save_lag$theta)/2
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
    mu <-  1/(1+exp(-X_O%*%res_save_post_lag$Thetas[k,1:p_O]))
    M  <- exp(X_OM%*%res_save_post_lag$Thetas[k,p_O+(1:p_OM)])
    Alpha_O[upper.tri(Alpha_O,diag=T)]  <- M*mu
    Beta_O[upper.tri(Alpha_O,diag=T)]   <- M*(1- mu)

    Reported_sample_post_lag <-fill.ReportBB.postlag_bi(res_save_post_lag$Death_est[k,],
                                            Alpha_O[1:pred_time,1:(pred_time +report.day)],
                                            Beta_O[1:pred_time, 1:(pred_time +report.day)],
                                            res_save_post_lag$pi_vec[k],
                                            res_save_post_lag$p_vec[k],
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
    data <- sim.data[     sim.data$state >= states_u[k]-10 -(14-d) &
                            sim.data$state < states_u[k] - (14-d) &
                            sim.data$days.left ==d &
                            is.na(sim.data$y) == F &
                            is.na(sim.data$yhat) == F,]
    cat('k:',k,' d:',d,' m:', median(data$y-data$yhat, na.rm=T),' mean:', mean(data$y-data$yhat, na.rm=T),'\n')
    index <- sim.data$state==states_u[k] &
      sim.data$days.left ==d
    if(is.na( median(data$y-data$yhat, na.rm=T))==F & d < 10)
      sim.data$yhat_mod[index] <- sim.data$yhat[index] + median(data$y-data$yhat, na.rm=T)
  }

}
write_fst(sim.data, file.path("data", "processed", "model_lag_refit.fst"))
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
write_fst(model_benchmark_bias, file.path("data", "processed", "model_benchmark_bias.fst"))

sim_today <- sim.data[state==max(states_u),]
alpha=0.1
sim_today<- sim_today[,.(quantile(yhat_mod,prob=c(1-alpha/2), na.rm=T),
                         quantile(yhat_mod,prob=c(alpha/2), na.rm=T),
                         quantile(yhat_mod,prob=c(0.5), na.rm=T)), by = .(date)]
names(sim_today) <- c("date","ci_upper","ci_lower","predicted_deaths")
write_fst(sim_today, file.path("data", "processed", "model_latest.fst"))


print(cbind(SCPRS[order(days.left,decreasing = T), .(mean(V2, na.rm=T),mean(V3, na.rm=T)), by = .(days.left)],benchmark[state>='2020-04-21', .(mean(SCRPS, na.rm=T)), by = .(days_left)]))
