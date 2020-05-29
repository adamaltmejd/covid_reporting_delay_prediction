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

start.predict.day=25
path.to.files <- file.path("data")
#files <- list.files('/Users/jonaswallin/Dropbox/temp/simulation_results_beta',
#                    pattern = "^param",
#                    full.names = TRUE)
files <- list.files(paste(path.to.files,"/simulation_results",sep=""),
                             pattern = "^param",
                             full.names = TRUE)
alpha_CI <- 0.1

index_files <- 1:(length(files)-1)
result <- readRDS(file.path("data", "processed", "processed_data.rds"))


cl <- parallel::makeCluster(nclust,setup_strategy = "sequential", outfile="")
doParallel::registerDoParallel(cl)

#for(i in index_files){
result_par <- foreach(i = index_files)  %dopar% {
  predicition_by_reported_day <-list()
  library(Matrix)
  library(tidyr)
  library(fst)
  library(stringr)
  load(files[i]) #res_save_lag

  cat('i=',i,", ", file=stdout())
  Date <- as.Date(result$dates_report[dim(res_save_lag$Death_est)[2]])
  cat('file= ',format(Date) , '\n',file=stdout())
  file.name = paste(path.to.files,"/simulation_results/prediction_valdiation.RData",sep="")
  write.data <- file.exists(file.name)==F
  if(write.data==F){
    load(paste(path.to.files,"/simulation_results/prediction_valdiation.RData",sep=""))
    write.data = Date%in%as.Date(names(result_par)) ==F
  }
  if(write.data){
    result <- readRDS(file.path("data", "processed", "processed_data.rds"))
    result$detected[14, 21:dim(result$detected)[2]] <- result$detected[14, 21:dim(result$detected)[2]] - 15

    pred_time = which(result$dates_report==str_extract(files[i],'[0-9]{4}-[0-9]{2}-[0-9]{2}'))
    cat('pred_time=',pred_time,"\n", file=stdout())
    N <- dim(result$detected)[1]
    lag <- res_save_lag$lag
    report <- splitlag(result$detected,as.Date(result$dates_report) ,res_save_lag$lag)

    X_T <- setup_data_lag(N, res_save_lag$npar, result$dates_report, res_save_lag$npar)
    X_O <- setup_data_postlag2(N, res_save_lag$lag, res_save_post_lag$npars, result$dates_report)
    Reported <- result$detected[1:pred_time,1:pred_time]
    Reported_fill <- cbind(Reported, matrix(NA, nrow=pred_time,ncol = (N-pred_time)))

    p_T <- length(res_save_lag$theta)/2
    p_O <- dim(X_O)[2]
    sim <- dim(res_save_lag$Death_est)[1]
    pred_set <-array(NA, dim = c(pred_time,N-pred_time,sim))

    Alpha_T <- matrix(NA, N,N)
    Beta_T  <- matrix(NA, N,N)

    mu <-  1/(1+exp(-X_T%*%res_save_lag$theta[1:p_T]))
    M  <- exp(X_T%*%res_save_lag$theta[p_T+(1:p_T)])
    Alpha_T[upper.tri(Alpha_T,diag=T)] <- M*mu
    Beta_T[upper.tri(Beta_T,diag=T)]   <- M*(1 - mu)

    Alpha_T <- Alpha_T[1:pred_time,1:N]
    Beta_T  <- Beta_T[1:pred_time, 1:N]
    Alpha_O <- matrix(NA, N,N)
    Beta_O  <- matrix(NA, N,N)
    mu <-  1/(1+exp(-X_O%*%res_save_post_lag$theta[1:p_O]))
    M  <- exp(X_O%*%res_save_post_lag$theta[p_O+(1:p_O)])
    Alpha_O[upper.tri(Alpha_O,diag=T)]  <- M*mu
    Beta_O[upper.tri(Alpha_O,diag=T)]   <- M*(1 - mu)

    for(k in 1:sim){


      Reported_sample_lag <-fill.ReportBB.lag(res_save_lag$Death_est[k,],
                                      Alpha_T,
                                      Beta_T,
                                      Reported_fill,
                                      result$dates_report[1:N],
                                      lag = lag)
      Reported_sample_post_lag <-fill.ReportBB.postlag(res_save_post_lag$Death_est[k,],
                                              Alpha_O,
                                              Beta_O,
                                              Reported_sample_lag,
                                              result$dates_report[1:N],
                                              lag = lag)
      pred_set[,,k] <- Reported_sample_post_lag[,(pred_time+1):N]
    }

    CI <- apply(pred_set,c(1,2), quantile,probs=c(alpha_CI/2,0.5,1-alpha_CI/2),na.rm=T)
    CI_low <- as.matrix(CI[1,,])
    colnames( CI_low)<- format(result$dates_report[(pred_time+1):N],"%Y-%m-%d")
    rownames( CI_low)<- format(result$dates[1:pred_time],"%Y-%m-%d")

    med <- as.matrix(CI[2,,])
    colnames( med)<- format(result$dates_report[(pred_time+1):N],"%Y-%m-%d")
    rownames( med)<- format(result$dates[1:pred_time],"%Y-%m-%d")

    CI_up <- as.matrix(CI[3,,])
    colnames( CI_up)<- format(result$dates_report[(pred_time+1):N],"%Y-%m-%d")
    rownames( CI_up)<- format(result$dates[1:pred_time],"%Y-%m-%d")
    index <- format(result$dates_report[pred_time],"%Y-%m-%d")

    predicition_by_reported_day$CI_low <- CI_low
    predicition_by_reported_day$CI_up <- CI_up
    predicition_by_reported_day$med <- med

    deaths_est_T <- apply(result$detected, 1, max, na.rm=T)
    #handling annoying misstakes

    predicition_by_reported_day$Truth <- as.matrix(result$detected[1:pred_time,(pred_time+1):N])
    colnames( predicition_by_reported_day$Truth)<- format(result$dates_report[(pred_time+1):N],"%Y-%m-%d")
    rownames( predicition_by_reported_day$Truth)<- format(result$dates[1:pred_time],"%Y-%m-%d")
    SCPRS <- matrix(0,nrow=pred_time, ncol = N-pred_time)
    for(k_death in 1:pred_time){
      if(k_death<=5)
        next
      for(k_report in 1:(N-pred_time)){

        predSample <-pred_set[k_death,k_report,]
        P1 <- predSample[1:ceiling(sim/3)]
        P2 <- predSample[(ceiling(sim/3) + 1):(floor(2*sim/3) )]
        P3 <- predSample[(floor(2*sim/3) + 1):(floor(2*sim/3) + length(P2))]
        EabsXX <- mean(abs(sample(P2)-sample(P3)))
        if(EabsXX == 0){
          EabsXX <- 1
        }
        Y <- predicition_by_reported_day$Truth[k_death,k_report]
        EabsYX <- mean(abs(Y-P1))
        SCPRS[k_death, k_report] <- -EabsYX/EabsXX - 0.5*log(EabsXX)
      }
    }
    colnames( SCPRS)<- format(result$dates_report[(pred_time+1):N],"%Y-%m-%d")
    rownames( SCPRS)<- format(result$dates[1:pred_time],"%Y-%m-%d")
    predicition_by_reported_day$SCPRS <- SCPRS
  }else{

    ttt<-which(names(result_par)==Date)
    cat('loading=',ttt,"\n", file=stdout())
    predicition_by_reported_day <- result_par[[ttt]]
  }
  predicition_by_reported_day
}

parallel::stopCluster(cl)
names(result_par) <- result$dates_report[25+index_files-1]
save(
     result_par,
     file=paste(path.to.files,"/simulation_results/prediction_valdiation.RData",sep="")
     )
nDelay <- 14
Pred_vec <- matrix(0,nrow=nDelay, ncol=7)
state=c()
date=c()
target = c()
days_left = c()
ci_upper = c()
ci_lower = c()
predicted_deaths = c()
SCRPS_out  = c()
for(day in names(result_par)){
  for(Delay in 1:(nDelay)){



    Report_day <- as.Date(day)+Delay
    Death_day <-  Report_day - nDelay
    if(Death_day<=as.Date('2020-04-06'))
      next

    CI_low <- result_par[[day]]$CI_low
    Row <- which(rownames(CI_low)==Death_day)
    Col <- which(colnames(CI_low)==Report_day)
    truth <- result_par[[day]]$Truth[Row,Col]
    if(length(Col)*length(Row)>0){
      CIl <- result_par[[day]]$CI_low[Row,Col]
      CIu <- result_par[[day]]$CI_up[Row,Col]

      med   <- result_par[[day]]$med[Row,Col]
      SCPRS   <- result_par[[day]]$SCPRS[Row,Col]
      state     = c(state,as.Date(day))
      date      = c(date, as.Date(Death_day))
      target    = c(target, truth)
      days_left = c(days_left,Report_day-as.Date(day))
      ci_upper  = c(ci_upper, CIu)
      ci_lower  = c(ci_lower, CIl)
      predicted_deaths = c(predicted_deaths,med)
      SCRPS_out  = c(SCRPS_out,SCPRS)
        Pred_vec[Delay,1] <- Pred_vec[Delay,1] + 1 #number
        Pred_vec[Delay,2] <- Pred_vec[Delay,2] + (truth >= CIl)*(truth <= CIu) #in CI
        Pred_vec[Delay,3] <- Pred_vec[Delay,3] + CIu - CIl #interval width
        Pred_vec[Delay,4] <- Pred_vec[Delay,4] + SCPRS #SCPRS
        Pred_vec[Delay,5] <- Pred_vec[Delay,5] + CIl - truth #in CI
        Pred_vec[Delay,6] <- Pred_vec[Delay,6] + CIu - truth  #in CI
        Pred_vec[Delay,7] <- Pred_vec[Delay,7] +  max(truth - CIu,CIl-truth,0)  #in CI

    }
  }
}

out_BGP <- data.table(state     = as.Date(state,origin="1970-01-01"),
                  date      = as.Date(date, origin="1970-01-01"),
                  target    = target,
                  days_left =  as.Date(state,origin="1970-01-01") - as.Date(date, origin="1970-01-01"),
                  ci_upper  = ci_upper,
                  ci_lower  = ci_lower,
                  predicted_deaths = predicted_deaths,
                  SCRPS   = SCRPS_out)

write_fst(out_BGP, file.path("data", "processed", "model_benchmark.fst"))
if(0){
  pdf('Prediction Result.pdf')
  par(mfrow=c(2,2))
  plot(1:(nDelay),Pred_vec[,2]/Pred_vec[,1],
       xlab='Days left to pred',
       ylab='90% CI coverage')
  plot(1:(nDelay),Pred_vec[,3]/Pred_vec[,1],
       xlab='Days left to pred',
       ylab='mean 90% CI width',
       type='l')
  plot(1:(nDelay),Pred_vec[,4]/Pred_vec[,1],
       xlab='Days left to pred',
       ylab='mean SCPRS',
       type='l')
  plot(1:(nDelay),Pred_vec[,7]/(Pred_vec[,1]-Pred_vec[,2]),
      xlab='Days left to pred',
      ylab='average error',
      type='l')
  dev.off()

  plot(1:(nDelay),Pred_vec[,5]/Pred_vec[,1],
       xlab='days left to obs',
       ylab='average centered CI',
       type='l',
       ylim=c(min(Pred_vec[,5]/Pred_vec[,1]), max(Pred_vec[,6]/Pred_vec[,1])))
  lines(1:(nDelay),Pred_vec[,6]/Pred_vec[,1])
  date_predict <- '2020-05-09'
  date_death   <- '2020-05-03'
  rbind(result_par[[date_predict]]$CI_low[date_death,],
        result_par[[date_predict]]$Truth[date_death,],
        result_par[[date_predict]]$CI_up[date_death,],
        result_par[[date_predict]]$SCPRS[date_death,])
  plot(result_par[[date_predict]]$SCPRS[date_death,])

  hist(pred_set[32,3,],breaks=60, main='2020-05-10 sann = 64',xlab='rapport',prob=T)
  hist(pred_set[32,3,],breaks=60, main='2020-05-10 sann = 64',xlab='rapport')

}
