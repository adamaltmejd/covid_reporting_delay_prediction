###
# benchmarking
#
###
library(foreach)
library(fst)
library(data.table)
source("src/util.r")
nclust <- 6
maxusage.day <- 20
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

Reported_T = result$detected
deaths_est_T <- apply(Reported_T, 1, max, na.rm=T)
data_T <- newDeaths(deaths_est_T,
                    Reported_T,
                    maxusage.day =maxusage.day)
cl <- parallel::makeCluster(nclust,setup_strategy = "sequential", outfile="")
doParallel::registerDoParallel(cl)

#for(i in index_files){
result_par <- foreach(i = index_files)  %dopar% {
  predicition_by_reported_day <-list()
  library( stringr)
  library(Matrix)
  library(tidyr)
  library(fst)

  load(files[i])
  cat('i=',i,", ", file=stdout())
  Date <- as.Date(str_extract(files[i],'[0-9]{4}-[0-9]{2}-[0-9]{2}'))
  cat('file= ',format(Date) , '\n',file=stdout())
  file.name = paste(path.to.files,"/simulation_results/prediction_valdiation.RData",sep="")
  write.data <- file.exists(file.name)==F
  if(write.data==F){
    load(paste(path.to.files,"/simulation_results/prediction_valdiation.RData",sep=""))
    write.data = Date%in%as.Date(names(result_par)) ==F
  }
  if(write.data){
    list2env(res_save, envir = environment())
    result <- readRDS(file.path("data", "processed", "processed_data.rds"))
    pred_time = which(result$dates_report==str_extract(files[i],'[0-9]{4}-[0-9]{2}-[0-9]{2}'))

    cat('pred_time=',pred_time,"\n", file=stdout())

    Reported_T = result$detected
    j <- pred_time
    j_0 <- max(1,j-7*3)
     Reported <- result$detected[j_0:j,j_0:j]
    dates_report <- result$dates_report[j_0:j]
    N_T <- dim(Reported_T)[1]
    deaths_est_T <- apply(Reported_T, 1, max, na.rm=T)
    data_T <- newDeaths(deaths_est_T,
                        Reported_T,
                        maxusage.day =maxusage.day)
    X_T <- setup_data(N_T, maxusage.day, result$dates_report, unique.days)

    #handling annoying misstakes
    data_mod <- newDeaths(deaths_est_T,
                        Reported_T,
                        maxusage.day =Inf)
    detected_mod<- result$detected
    for(k in 1:dim(data_T$report.new)[1])
      detected_mod[k,] <- cumsum(replace_na(data_mod$report.new[k,],0))

    detected_mod[is.na(result$detected)] <- NA


    Reported <- result$detected[j_0:pred_time,j_0:pred_time]
    Reported_fill <- cbind(Reported, matrix(NA, nrow=(pred_time-j_0+1),ncol = (N_T-pred_time)))

    p <- dim(Thetas)[2]/2
    sim <- dim(Thetas)[1]
    pred_set <-array(NA, dim = c(pred_time-j_0+1,N_T-pred_time,sim))
    for(k in 1:sim){
      beta_1 <- Thetas[k,1:p]
      beta_2 <- Thetas[k,(p+1):(2*p)]
      Alpha_T <- matrix(NA, N_T,N_T)
      Beta_T  <- matrix(NA, N_T,N_T)

      Alpha_T[upper.tri(Alpha_T,diag=T)] <- exp(X_T%*%beta_1)
      Beta_T[upper.tri(Beta_T,diag=T)]   <- exp(X_T%*%beta_2)
      Alpha_T <- Alpha_T[j_0:j,j_0:N_T]
      Beta_T  <- Beta_T[j_0:j, j_0:N_T]
      Reported_sample <-fill.ReportBB(Death_est[k,],
                                      Alpha_T,
                                      Beta_T,
                                      Reported_fill,
                                      maxusage.day = true.day)
      pred_set[,,k] <- Reported_sample[,(pred_time+2-j_0):(N_T-j_0+1)]
    }

    CI <- apply(pred_set,c(1,2), quantile,probs=c(alpha_CI/2,0.5,1-alpha_CI/2))
    CI_low <- as.matrix(CI[1,,])
    colnames( CI_low)<- format(result$dates_report[(pred_time+1):N_T],"%Y-%m-%d")
    rownames( CI_low)<- format(result$dates[j_0:pred_time],"%Y-%m-%d")

    med <- as.matrix(CI[2,,])
    colnames( med)<- format(result$dates_report[(pred_time+1):N_T],"%Y-%m-%d")
    rownames( med)<- format(result$dates[j_0:pred_time],"%Y-%m-%d")

    CI_up <- as.matrix(CI[3,,])
    colnames( CI_up)<- format(result$dates_report[(pred_time+1):N_T],"%Y-%m-%d")
    rownames( CI_up)<- format(result$dates[j_0:pred_time],"%Y-%m-%d")
    index <- format(result$dates_report[pred_time],"%Y-%m-%d")

    predicition_by_reported_day$CI_low <- CI_low
    predicition_by_reported_day$CI_up <- CI_up
    predicition_by_reported_day$med <- med
    predicition_by_reported_day$Truth <- as.matrix(detected_mod[j_0:pred_time,(pred_time+1):N_T])

    colnames( predicition_by_reported_day$Truth)<- format(result$dates_report[(pred_time+1):N_T],"%Y-%m-%d")
    rownames( predicition_by_reported_day$Truth)<- format(result$dates[j_0:pred_time],"%Y-%m-%d")
    SCPRS <- matrix(nrow=pred_time-j_0+1, ncol = N_T-pred_time)
    for(k_death in 1:(pred_time-j_0+1)){
      for(k_report in 1:(N_T-pred_time)){

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
    colnames( SCPRS)<- format(result$dates_report[(pred_time+1):N_T],"%Y-%m-%d")
    rownames( SCPRS)<- format(result$dates[j_0:pred_time],"%Y-%m-%d")
    predicition_by_reported_day$SCPRS <- SCPRS
  }else{

    ttt<-which(names(result_par)==Date)
    cat('loading=',ttt,"\n", file=stdout())
    predicition_by_reported_day <- result_par[[ttt]]
  }
  predicition_by_reported_day
}

parallel::stopCluster(cl)
start.predict.day=40
names(result_par) <- result$dates_report[start.predict.day+index_files-1]
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
