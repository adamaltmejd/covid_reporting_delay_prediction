####
# data wrangling
#
#####
library(Matrix)
library(data.table)
library(rSPDE)
#'  returns a vector of how many has deaths are reported at t days after
#'  @param t       - (int > 0 ) how many days
#'  @param reports - (N x N) reported cumlative deaths
deaths_at_t <-function(reports, t){

    N <- dim(reports)[1]
    dt <- diag(reports[,-(1:t)])
    dt <- c(dt, rep(NA, N- length(dt)))
    names(dt) <- rownames(reports)
    return(dt)
}




##
# how m
# transforms data,
#' we remove data if negative..
#'
#'  @param deaths  - (N x 1) how many has died (thruth)
#'  @param reports - (N x N) reported cumlative deaths
#'  @param maxusage.day - (int) only use data up to maxusage.days i.e.
#'                       reports[i,i + maxusage.day - 1]
#' @return list
#'             - death.remain how many remains to be reported
#'             - report.new   how many are reported
##
newDeaths <-function(deaths, reports,maxusage.day = -1){
    newreport <- reports
    #reports_temp[is.na(reports_temp)]=0
    ##
    # ugly fix
    ##

    reports_temp <- reports

    for(i in 1:dim(reports_temp)[1]){
        reports_temp[i,1:(i-1)]=0
        if(is.na(reports_temp[i,i]))
            reports_temp[i,i] <- 0
        for(j in i:dim(reports_temp)[2]){
            if(is.na(reports_temp[i,j]))
                reports_temp[i,j]= reports_temp[i,j-1]
        }

    }
    reports_temp <- cbind(0,reports_temp)
    newreport <- t(diff(t(reports_temp)))
    newreport[is.na(reports)]=NA
    #newreport[upper.tri(newreport)] <- diff.report[upper.tri(diff.report,T)]
    newreport[newreport<0 & is.na(newreport)==F]=0 #fake
    death.rem <- diag(deaths)

    cumsum_report <- t(apply(newreport, 1, function(x){x[is.na(x)==F]=cumsum(x[is.na(x)==F]); return(x)}))
    dr<-deaths - cumsum_report
    N <- length(deaths)
    dr <- dr[1:(N-1),1:(N-1)]
    death.rem[upper.tri(newreport)] <- dr[upper.tri(dr,T)]
    death.rem[lower.tri(death.rem)] <- NA
    diag(death.rem)[is.na(diag(reports))] <- NA

    ##
    # ugly fix 2
    ##
    for(k in 1:dim(death.rem)[1]){
        index <- is.na(death.rem[k,]) & is.na(newreport[k,])==F
        if(sum(index)>0){
            death.rem[k, index] <- deaths[k]
        }
    }
    if(maxusage.day>0){
        for(i in 1:length(death.rem)){
            if(i + maxusage.day - 1 < N ){
                death.rem[i,(i+maxusage.day):N] = NA
                newreport[i,(i+maxusage.day):N] = NA
            }
        }
    }
    #removing small bugs in reporting
    newreport[death.rem <0 & is.na(death.rem)==F] = 0
    death.rem[death.rem <0& is.na(death.rem)==F] = 0
    index <- (death.rem< newreport) &  is.na(death.rem) ==F & is.na(newreport)==F
    newreport[index] =death.rem[index]
    return(list(death.remain = death.rem, report.new = newreport))
}


##
# data exploration build data to explore
#
#' @param  (list) $report.new , $death.remain
#' @parma  (m x 1) how many lags should be collected
##
data.by.lag <- function(Data, lags = 1:5){


    d <- dim(Data$report.new)[1]
    lag <- 1
    result <- c()
    TotalD <- diag(Data$death.remain)
    for(lag in lags){
        n = diag(Data$death.remain[,-(1:lag)])
        N = TotalD[1:(d-lag)]
        y =  diag(Data$report.new[,-(1:lag)])
        ys <- y + (N-n)
        date.reported <- colnames(Data$report.new)[(lag+1):d]
        date          <- colnames(Data$report.new)[1:(d-lag)]
        index = is.na(y)==F
        if(sum(index)>0){
            dat.temp = data.table(
                                    date           = as.Date(date[index]),
                                    date.reported  = as.Date(date.reported[index]),
                                    lag            = lag,
                                    n              = n[index],
                                    y              = y[index],
                                    ys             = ys[index],
                                    N              = N[index])
            result <- rbind(result,dat.temp)
        }
    }
    return(result)
}

####
# model functions
#
####




##
# build holiday covariates vector for holiday
#
##
#  holidays - (N x 1) true if day is holiday false else
##
buildXholiday <- function(N,holidays){
    ##
    # base matrix
    ##
    index_base <- t(matrix(rep(holidays,N),ncol = N))

    index_base <- index_base[upper.tri(index_base,diag=T)]
    index_base <- 1* index_base
    #sparse matrix
    i_base <- 1:length(index_base)
    i_base <- i_base[index_base==1]
    return(sparseMatrix(j=rep(1,length(i_base)),i=i_base,  dims=c(length(index_base), 1)))
}


##
# build day effect matrix
#  nDayEffects - number of speical days effect (1- first day, 2- first + second day)
#  N - number of days
#  nDayEffects - how many day covariate effect to create
##
buildXdayeffect <- function(N, nDayEffects = 1){
    nDayEffects <- min(N,nDayEffects)
    index_days <- toeplitz(c( (1:nDayEffects), rep(0,N -nDayEffects)))
    index_days <- index_days[upper.tri(index_days,diag=T)]
    j_ <-  index_days[index_days>0]
    i_base <- 1:length(index_days)
    i_ <-  i_base[index_days>0]
    return(sparseMatrix(i=i_,j=j_ , dims=c(length(index_days), nDayEffects)  ) )
}

##
#'
#' builds the corresponding X matrix for the upper triangular matrix
#'
#' @param report (n x n) reported day 1 not 0
#' @parma nrepporteffect  (int) how long lag
#'
buildXreporteffect <- function(report, nreporteffect = 1){

    n = dim(report)[2]
    reportNotNa <- report
    reportNotNa[is.na(reportNotNa)] = 0
    X_base <- t(apply(reportNotNa,1, cumsum))
    X_base[is.na(report) | report==0] = -1

    X_base <- X_base[upper.tri(X_base,diag=T)]
    j_ <- X_base[X_base>=0 & X_base < nreporteffect]
    j_ <- j_ + 1 # zero is included
    i_base <- 1:length(X_base) # position on the upper triangular matrix
    i_ <-  i_base[X_base>=0 & X_base < nreporteffect]

    Res <- sparseMatrix(i=i_,j=j_ , dims=c(length(X_base), nreporteffect )  )
    return(Res)
}

buildXNoReportLagOne <- function(report){

    ReportLagOne = c(diag(report[,-c(1)]),0)==0
    X_base <- matrix(rep(ReportLagOne, length(report[,1])), nrow=length(report[,1]) )
    X_base <- X_base[upper.tri(X_base,diag=T)]
    return(X_base)
}
buildXNoreportPrev <- function(report){

    n = dim(report)[2]
    reportNotNa <- report
    reportNotNa[is.na(reportNotNa)] = 1
    X_base <- t(apply(reportNotNa==0,1, function(x) {c(0,cumsum(x))[1:length(x)]}))

    X_base[cbind(1,reportNotNa[,1:(n-1)])==1] = 0

    X_base <- X_base[upper.tri(X_base,diag=T)]
    i_ <- 1:length(X_base) # position on the upper triangular matrix
    i_ <-  i_[X_base>0]
    j_ <- rep(1,length(i_))

    Res <- sparseMatrix(i=i_,j=j_, x = X_base[X_base>0] , dims=c(length(X_base), 1 )  )
    return(Res)
}


###
# setting up the fixed effect matrix
#
#
##
setup_data.uk <- function(dates_report, Predict.day, unique.prob=NULL, reported =NULL){
    N <- length(dates_report)
    if(is.null(reported)){
        X <- buildXdayeffect(N,Predict.day)
    }else{
        X <- buildXreporteffect(reported, Predict.day)
    }
    #holidays and weekends

    #put together non unique days
    if(is.null(unique.prob)==F)
        X <- cbind(X[,1:unique.prob,drop=F],rowSums(X[,1:unique.prob,drop=F])==0)

    unique.weekdays <- unique( weekdays(dates_report))
    X.days <- c()
    for(i in 1:length(unique.weekdays)){
        day.index <- weekdays(dates_report)%in%unique.weekdays[i]
        X.day <-  buildXholiday(N, day.index)
        X.days <- cbind(X.days, X.day)
    }
    index.shift <- dates_report>=as.Date("2021-01-27")
    X.shift <-  buildXholiday(N, index.shift)
    X_dim <- dim(X)[2]
    X <- cbind(X, X.days, X.shift)
    colnames(X) <- c(paste('lag_',0:(X_dim-1),sep=''), unique.weekdays, 'after 2021-01-27')
    return(X)


}
###
# setting up the fixed effect matrix
# for beta
##
X.uk <- function(dates_report, reported){
    X <- setup_data.uk(dates_report, length(dates_report), 5)
    LagOneNoReport <- buildXNoReportLagOne(reported)
    X_lag <- X[,14] * LagOneNoReport
    X[,14] <- X_lag*X[,3]
    colnames(X)[14] <- "lag_2 after 2021-01-27"
    Sunday <- X[,12]
    X[,12] <- X[,13]+X[,12]
    colnames(X)[12] <- "monday or sunday"
    X[,13] <- X[,13]*X[,4]
    colnames(X)[13] <- "monday lag 4"
    X <- cbind(X, (X[,3] + X[,2])*X[,12] )
    colnames(X)[15] <- "monday or sunday lag 1 or 2"
    X <- cbind(X, X_lag)
    colnames(X)[16] <- "after 2021-01-27"
    X <- cbind(X,X[,14]*X[,3]*Sunday)
    colnames(X)[17] <- "sunday lag_2"
    X <- X[,c(2,3,4,12,13,14,15,16,17)]

    X_temp <- cbind(1,X)
    colnames(X_temp) <- c("intercept",colnames(X))
    X <- X_temp
    X_rep <- buildXNoreportPrev(reported)
    X_temp <- cbind(X, X_rep)
    colnames(X_temp) <- c(colnames(X),"no rep")
    X <- X_temp
    return(X)


}

###
# setting up the fixed effect matrix
# for beta
##
X.uk.reported <- function(dates_report, reported){
    X <- setup_data.uk(dates_report, length(dates_report), 5, reported)

    X <- X[,c(2,3,4,12,13)]
    X <- cbind(1,X)
    colnames(X)[1] <- "intercept"

    X_rep <- buildXNoreportPrev(reported)
    X <- cbind(X, X_rep)
    colnames(X)[dim(X)[2]] <- "no prev rep"
    return(X)


}

#'
#' regresion for zero inflation
#'
#'
X.uk.zero.reported <- function(dates_report, reported){
    X <- setup_data.uk(dates_report, max.days.to.report, 5, reported)
    X_lag <- X[,14]
    X[,14] <- X[,14]*X[,3]
    colnames(X)[14] <- "lag_3 after 2021-01-27"
    X[,12] <- X[,13]+X[,12]
    colnames(X)[12] <- "monday or sunday"
    X[,13] <- X[,13]*X[,4]
    colnames(X)[13] <- "monday lag 4"
    X <- cbind(X, (X[,3] + X[,2])*X[,12] )
    colnames(X)[15] <- "monday or sunday lag 2 or 3"
    X <- cbind(X, X_lag)
    colnames(X)[16] <- "after 2021-01-27"
    X <- X[,c(2,3,4,12,13,14,15,16)]


    X <- cbind(1,X)
    colnames(X)[1] <- "intercept"
    X <- X[,1 , drop=F]
    return(X)


}
#'
#' regresion for zero inflation
#'
#'
X.uk.zero <- function(dates_report){
    X <- setup_data.uk(dates_report, max.days.to.report, 5)
    X_lag <- X[,14]
    X[,14] <- X[,14]*X[,3]
    colnames(X)[14] <- "lag_3 after 2021-01-27"
    X[,12] <- X[,13]+X[,12]
    colnames(X)[12] <- "monday or sunday"
    X[,13] <- X[,13]*X[,4]
    colnames(X)[13] <- "monday lag 4"
    X <- cbind(X, (X[,3] + X[,2])*X[,12] )
    colnames(X)[15] <- "monday or sunday lag 2 or 3"
    X <- cbind(X, X_lag)
    colnames(X)[16] <- "after 2021-01-27"
    X <- X[,c(2,3,4,12,13,14,15,16)]


    X <- cbind(1,X)
    colnames(X)[1] <- "intercept"
    X <- X[,1 , drop=F]
    return(X)


}



###
# optimizing to find coeffients for probabilites of reporting
# i.e. mu, and M in beta binomial distribution
#
#' @param result - list
#'                      $detected (N x N) number of deatected detected dead
#'                      $report   (N x N) was date reported at report.date
#'                      $dates    (N x 1) which Date
#' @param max.days.to.report -  after max.days.to.report data is assume known
#' @param zero.inflation  - should the model be zero infaltion
###
fit.mu.M <- function(result, max.days.to.report, zero.inflation=FALSE, use.reports = T){


    N <- deaths_at_t(result$detected, max.days.to.report)
    fixed.dates <- sum(is.na(N)==F)
    N.fixed <- N[1:fixed.dates]

    Data <- newDeaths(N.fixed,
                      reports = result$detected[1:fixed.dates,1:fixed.dates ],
                      maxusage.day = max.days.to.report)

    dates_report <- result$dates_report[1:fixed.dates]
    reports <- result$report[1:fixed.dates,1:fixed.dates]
    if(use.reports){
        X_T <- X.uk.reported(dates_report,reports)
    }else{
        X_T <- X.uk(dates_report,reports)
    }
    if(zero.inflation){
        if(use.reports){
            X_zero <- X.uk.zero.reported(dates_report,reports)
        }else{
            X_zero <- X.uk.zero(dates_report)
        }
    }
    index <- upper.tri(Data$death.remain,diag = T)
    Data$report.new[ is.na(reports)==F &reports==0 ] = NA
    y = Data$report.new[index]
    n = Data$death.remain[index]
    index = is.na(y)==F
    X <- as.matrix(X_T[index,])
    beta_index <- colSums(X)>10
    X <- X[,beta_index]
    if(zero.inflation)
        X_zero <- as.matrix(X_zero[index,])
    y <- y[index]
    n <- n[index]
    #inital fit

    fitShort <- glm( cbind(y, n) ~ -1 + X,
                     family = "binomial")
    theta0.mu  <- c(fitShort$coefficients)
    theta0.M   <- rep(0, length(theta0.mu))
    X1 <- X
    X2 <- X
    if(zero.inflation){
        X3 <- X_zero
        theta0.pi <- rep(0, dim(X_zero)[2])
        lik <- function(x){ - log_bb.mixed.zero(x, y, n, X1, X2,X3)}
        res.optim <- optim(c(theta0.mu ,theta0.M, theta0.pi), lik)
        #fixe
        res.optim <- optim(res.optim$par, lik)
        res.optim <- optim(res.optim$par, lik)
        res.optim <- optim(res.optim$par, lik)
        res.optim <- optim(res.optim$par, lik)
        res.optim <- optim(res.optim$par, lik)
        res.optim <- optim(res.optim$par, lik)
        p1 <- dim(X1)[2]
        p2 <- dim(X2)[2]
        p3 <- dim(X3)[2]
        beta <- list(mu1 = res.optim$par[1:p1] ,
                     M1 = res.optim$par[(p1+1):(p1 + p2)],
                     pi  =res.optim$par[(p1+p2+1):(p1 + p2 + p3)] )
        names(beta$M1) <- names(beta$mu1)
        names(beta$pi) <- colnames(X3)
    }else{
        lik <- function(x){ - log_bb(x, y, n, X1, X2)}
        res.optim <- optim(c(theta0.mu,theta0.M), lik)
        res.optim <- optim(res.optim, lik)
        res.optim <- optim(res.optim, lik)
        beta <- list(mu = res.optim$par[1:dim(X1)[2]] ,
                     M  =res.optim$par[(dim(X1)[2]+1):(dim(X1)[2] + dim(X2)[2])] )
        names(beta$M) <- names(beta$mu)
    }
    beta$index = beta_index
    return(beta)
}

#'
#' @param result - list
#'                      $detected (N x N) number of deatected detected dead
#'                      $report   (N x N) was date reported at report.date
#'                      $dates    (N x 1) which Date
#' @param max.days.to.report -  after max.days.to.report data is assume known
#' @param samples            - (int) how many mcmcmc samples
#' @param burnin.perc        - ([0,1]) how many precetentage of samples should be burning
#' @oaram samples.local      - (int) samples in uniform HMR within MCCMC
sample.uk.deaths <- function(result,
                             max.days.to.report,
                             samples,
                             burnin.perc = 0.1,
                             samples.local=10 ,
                             zero.inflation  = TRUE){


    Nest <- apply(result$detected,1, max, na.rm=T)
    data_full <- newDeaths(Nest, reports = result$detected, maxusage.day = max.days.to.report)
    n.days <- length(result$dates)
    beta <- fit.mu.M(result, max.days.to.report, zero.inflation, use.reports = F)
    X_full <- X.uk(result$dates_report, result$report)
    X_full <- X_full[,beta$index] #remove emptty covariates
    p <- dim(X_full)[2]
    mu <- 1/(1+exp(-X_full%*%beta$mu1 ))
    M  <- exp(X_full%*%beta$M1)

    Alpha <- matrix(NA,
                    ncol=n.days,
                    nrow= n.days)
    Beta  <- Alpha
    if(zero.inflation){
        X_zero <- X.uk.zero(result$dates_report)
        Prob   <- Alpha
        Prob[upper.tri(data_full$report.new,diag=T)] <- 1/(1+exp(-X_zero%*%beta$pi))
    }

    Alpha[upper.tri(data_full$report.new,diag=T)] <- M * mu
    Beta[upper.tri(data_full$report.new,diag=T)]  <- M * (1-mu)


    data_full$report.new[ is.na(result$report)==F &result$report==0 ] = NA

    alpha.MCMC <- rep(ceiling(0.1*Nest), n.days)
    n.burnin <- ceiling(burnin.perc* samples)
    #burnin
    for(i in 1:(n.burnin)){
        res <- sample.deathsBB(Nest,
                               Alpha,
                               Beta,
                               data_full$report.new,
                               samples=samples.local,
                               alpha.MCMC,
                               max.days.to.report = max.days.to.report,
                               prob = Prob)
        alpha.MCMC[res$acc/samples.local > 0.3] <- ceiling(alpha.MCMC[res$acc/samples.local > 0.3] * 1.1)
        alpha.MCMC[res$acc/samples.local < 0.3] <- ceiling(alpha.MCMC[res$acc/samples.local < 0.3] * 0.9)
        Nest <- res$deaths
    }
    #samples
    N.samples <- matrix(NA, nrow= samples, ncol = n.days)
    colnames(N.samples) <- as.character(result$dates)
    result$detected[ is.na(result$report)==F &result$report==0 ] = NA
    for(i in 1:(samples)){
        res <- sample.deathsBB(Nest,
                               Alpha,
                               Beta,
                               data_full$report.new,
                               samples=samples.local,
                               alpha.MCMC,
                               max.days.to.report = max.days.to.report,
                               prob = Prob)
        N.samples[i,] <- res$deaths
    }
    return(N.samples)
}

###
# beta binomial components
#
###

#'
#' log likelihood function for beta binomial
#'
#' @param  theta (p1+p2 x 1 ) mu, M parameters for beta binomial
#' @param  y (k x 1)        number of positve
#' @param  n (k x 1)        number of observations
#' @param X1 (k x p1)      fixed effects fo Mu
#' @param X2 (k x p2)      fixed effects fo Mu
log_bb<- function(theta, y, n, X1, X2){
    p1 <- dim(X1)[2]
    p2 <- dim(X2)[2]
    theta_mu  <- theta[1:p1]
    theta_M   <- theta[(p1+1):(p1+p2)]
    mu <- 1/(1+exp(-X1%*%theta_mu))
    M  <- exp(X2%*%theta_M)
    alpha <- mu*M
    beta  <- M-alpha
    if(min(1/(alpha+beta)) <10^-3)
        return(-Inf)
    lik <- sum(dBB(x =y,
                   size = n,
                   alpha = alpha,
                   beta = beta,
                   log.p=T))
    return(lik)
}




#'
#' log likelihood function for mixture of two beta binomial bimodal
#'
#' @param  theta (2*(p1+p2) + p3 x 1 ) mu, M parameters for beta binomial
#' @param  y (k x 1)        number of positve
#' @param  n (k x 1)        number of observations
#' @param X1 (k x p1)      fixed effects of Mu
#' @param X2 (k x p2)      fixed effects of M
#' @param X3 (k x p3)      fixed effects of pi
log_bb.mixed<- function(theta, y, n, X1, X2, X3){
    p1 <- dim(X1)[2]
    p2 <- dim(X2)[2]
    p3 <- dim(X3)[2]
    theta_mu1  <- theta[1:p1]
    theta_mu2  <- theta[(p1+1):(2*p1)]
    theta_M1   <- theta[(2*p1+1):(2*p1 + p2)]
    theta_M2   <- theta[(2*p1+p2+1):(2*p1 + 2*p2)]
    theta_pi   <- theta[(2*p1+2*p2+1):(2*p1 + 2*p2 + p3)]
    mu1 <- 1/(1+exp(-X1%*%theta_mu1))
    M1  <- exp(X2%*%theta_M1)
    alpha1 <- mu1*M1
    beta1  <- M1-alpha1


    if(min(1/(alpha1+beta1)) <10^-3)
        return(-Inf)
    mu2 <- 1/(1+exp(-X1%*%theta_mu2))
    M2  <- exp(X2%*%theta_M2)
    alpha2 <- mu2*M2
    beta2  <- M2-alpha2

    if(min(1/(alpha2+beta2)) <10^-3)
        return(-Inf)

    pi =1/(1+exp(-X3%*%theta_pi))
    lik <- pi * dBB(x =y,
                    size = n,
                    alpha = alpha1,
                    beta = beta1,
                    log.p=F) +
           (1 -pi) * dBB(x =y,
                         size = n,
                         alpha = alpha2,
                         beta = beta2,
                         log.p=F)
    lik.sum <- sum(log(lik))
    return(lik.sum)
}




#'
#' log likelihood function beta binomial bimodal zero inflated
#' @param  theta (2*(p1+p2) + p3 x 1 ) mu, M parameters for beta binomial
#' @param  y (k x 1)        number of positve
#' @param  n (k x 1)        number of observations
#' @param X1 (k x p1)      fixed effects of Mu
#' @param X2 (k x p2)      fixed effects of M
#' @param X3 (k x p3)      fixed effects of pi
log_bb.mixed.zero <- function(theta, y, n, X1, X2, X3){
    p1 <- dim(X1)[2]
    p2 <- dim(X2)[2]
    p3 <- dim(X3)[2]
    theta_mu1  <- theta[1:p1]
    theta_M1   <- theta[(p1+1):(p1 + p2)]
    theta_pi   <- theta[(p1+p2+1):(p1 + p2 + p3)]
    mu1 <- 1/(1+exp(-X1%*%theta_mu1))
    M1  <- exp(X2%*%theta_M1)
    alpha1 <- mu1*M1
    beta1  <- M1-alpha1


    if(min(1/(alpha1+beta1)) <10^-3)
        return(-Inf)


    pi =1/(1+exp(-X3%*%theta_pi))
    lik <- pi * dBB(x =y,
                    size = n,
                    alpha = alpha1,
                    beta = beta1,
                    log.p=F) +
        (1 -pi) * (1* (y==0))

    lik.sum <- sum(log(lik))
    # regularization
    lik.sum <- lik.sum + sum(-theta_M1^2/2)
    return(lik.sum)
}

##
# explore if the distribution fits the data
#
#' @param  (list) $report.new , $death.remain
#' @param  Alpha [mxm] dbbeta - coeff
#' @param  Beta  [mxm] dbbeta - coeff
#' @param  Pi    [mxx] mixing coeff
#' @parmam lag (m x 1) how many lags should be collected
##
zero.BB.dist.by.lag <- function(Data,
                                Alpha,
                                Beta,
                                Pi,
                                lags = 1:5,
                                sim=1000,
                                probs = c(0.05,0.95)){


    d <- dim(Data$report.new)[1]
    lag <- 1
    result <- c()
    for(lag in lags){
        n = diag(Data$death.remain[,-(1:lag)])
        y =  diag(Data$report.new[,-(1:lag)])
        index = is.na(y)==F
        if(sum(index) > 0){
            alpha <- diag(Alpha[,-(1:lag)])
            beta  <- diag(Beta[,-(1:lag)])
            pi    <- diag(Pi[,-(1:lag)])
            size.y <- length(y)
            y.sample <- matrix(NA, ncol = size.y, nrow = sim)

            for(i in 1:sim){
                U <- runif(size.y) < pi
                p <- rbeta(size.y, alpha, beta)
                y <- rbinom(size.y, n, p)
                y[U==F] = 0
                y.sample[i, ] <- y
            }
            y.quant <- apply(y.sample,2, function(x) {quantile(x, probs=probs)})
            p0 = apply(y.sample,2, function(x) {mean(x==0)})
            date.reported <- colnames(Data$report.new)[(lag+1):d]
            date          <- colnames(Data$report.new)[1:(d-lag)]

            dat.temp = data.table(
                date = as.Date(date[index]),
                date.reported = as.Date(date.reported[index]),
                lag = lag,
                prob0 = p0[index],
                n     = n[index],
                low = y.quant[1,index],
                upp = y.quant[2,index])
            result <- rbind(result,dat.temp)
        }
    }
    return(result)
}



##
# PRobability of beta binomial
##
dBB<- function(x, size, alpha, beta, log.p = F){
    const = lgamma(size + 1) - lgamma(x + 1) -
        lgamma(size - x + 1)
    ld <- const +
        lgamma(x + alpha) + lgamma(size - x + beta) -
        lgamma(size + alpha + beta) +
        lgamma(alpha + beta) -
        lgamma(alpha) - lgamma(beta)
    if(log.p==F){
        return(exp(ld))
    }
    return(ld)
}


##
# sample reported deaths (N) after rep.day (if Inf) predicit day using Beta binomial dist
#
#'  @param  deaths      - (N x 1) estimate of number of deaths
#'  @param  alpha       - (N x N) matrix of beta binom parameter (only upper triangular part relevant)
#'  @param  beta        - (N x N) matrix of beta binom parameter (only upper triangular part relevant)
#'  @param  Reported    - (N x N) matrix of reported deaths cumlative (only upper triangular relevant)
#'  @param  samples     - (int) how many samples should be done for each fixed data  (cheap)
#'  @param  alpha.MCMC  - (N x 1) stepsizes for the discrete MH-RW
#' @param   max.days.to.report -  after max.days.to.report data is assume known
#' @param  prob         - (N x N) probability of zero inflation
##
sample.deathsBB <- function(deaths,
                            alpha,
                            beta,
                            Reported,
                            samples,
                            alpha.MCMC,
                            max.days.to.report,
                            prob = NULL){

    N <- length(deaths)
    acc <- rep(0,N)
    for(i in 1:N){
        if(i + max.days.to.report >= N){
            alpha_i     = alpha[i,i:N]
            beta_i      = beta[i,i:N]
            Reported_i  = Reported[i,i:N]
            index = is.na(Reported_i)==F
            alpha_i  = alpha_i[index]
            beta_i  = beta_i[index]
            if(is.null(prob)==F){
                prob_i     = prob[i,i:N]
                prob_i = prob_i[index]
            }else{
                prob_i = NULL
            }
            Reported_i  = Reported_i[index]
            if(length(Reported_i)>0){
                lik_i <- loglikDeathsGivenProbBB(deaths[i],alpha_i , beta_i, Reported_i, prob_i)
            }else{
                deaths[i] <- NA
                acc[i]    <- 0
                next
            }

            for(j in 1:samples){
                death_star <- sample((deaths[i]-alpha.MCMC[i]):(deaths[i]+alpha.MCMC[i]), 1)
                lik_star <- loglikDeathsGivenProbBB(death_star, alpha_i , beta_i, Reported_i, prob_i)
                if(is.nan(lik_star))
                    next
                if(log(runif(1)) < lik_star-lik_i){
                    lik_i = lik_star
                    deaths[i] <- death_star
                    acc[i] = acc[i] + 1
                }
            }
        }else{
            deaths[i] <- deaths[i]
        }
    }
    return(list(deaths=deaths, acc = acc))
}

##
# log liklihood of obseving report given death and prob
# density is Beta binomial
#  deaths  - (int) true number of deaths
#  alpha       - (n x 1) bb parameter 1
#  beta        - (n x 1) bb parameter 2
#  prob        - (n x 1) probability of not zero
#  report  - (n x 1) reported deaths each day
##
loglikDeathsGivenProbBB <- function(death, alpha, beta, report, prob = NULL){

    if(death < sum(report,na.rm=T))
        return(-Inf)
    n <- length(alpha)
    if(n>1){
        remain = death-cumsum(c(0, report[1:(n-1)]))
    }else{
        remain = death
    }

    if(is.null(prob)){
    return(sum(dBB(x =report,
                   size = remain,
                   alpha = alpha[is.na(alpha)==F],
                   beta = beta[is.na(alpha)==F],
                   log.p=T)))
    }else{
        dens <- prob[is.na(alpha)==F] * dBB(x =report,
                                            size = remain,
                                            alpha = alpha,
                                            beta = beta,
                                            log.p=F) +
                (1- prob[is.na(alpha)==F]) * (1* (report==0))
        return(sum(log(dens)))
    }
}
##
# Generates prediction data data.table with
# N mean median 0.05, 0.95 quantiles
#
#' @param result - list
#'                      $detected (N x N) number of deatected detected dead
#'                      $report   (N x N) was date reported at report.date
#'                      $dates    (N x 1) which Date
#'
#' @param   max.days.to.report -  after max.days.to.report data is assume known
#' @param report.dates  - which dates should be data generated for (if null) just latest date
##
uk.prediction <- function(result, max.days.to.report, report.dates =NULL){


    if(is.null(report.dates)){
        report.dates = max(result$dates_report)
    }
    deaths <- data.table(
                        date.reported = structure(numeric(0), class = "Date"),
                        date          = structure(numeric(0), class = "Date"),
                        mean          = numeric(),
                        median        =  numeric(),
                        low           =  numeric(),
                        upp           =  numeric())

    for(i in 1:length(report.dates)){
        result_i  <- result
        j_i <-  which(result_i$dates_report == report.dates[i])
        result_i$dates_report       <- result_i$dates_report[1:j_i]
        result_i$dates              <- result_i$dates[1:j_i]
        result_i$dates_not_reported <- result_i$dates_not_reported[1:j_i]
        result_i$report             <- result_i$report[1:j_i, 1:j_i]
        result_i$detected           <- result_i$detected[1:j_i, 1:j_i]

        N.est <- sample.uk.deaths(result_i, max.days.to.report, samples = 1000)
        N.quantile <- t(apply(N.est, 2, function(x) {c(mean(x, na.rm=T), quantile(x,probs=c(0.05,0.5,0.95), na.rm=T))}))
        deaths.i <- data.table(
            date.reported = report.dates[i],
            date          = as.Date(colnames(N.est)),
            mean          =  N.quantile[,1],
            median        =  N.quantile[,3],
            low           =  N.quantile[,2],
            upp           =  N.quantile[,4]
        )
        deaths <- rbind(deaths, deaths.i)
    }
    return(deaths)
}



####
# GP stuff
###

###
#' @param x parameter matern log(kappa,nu,sigma,sigma_noise)
#' @param D distance matrix
#' @param y data
##
likelihood_death <- function(x, D, y){
    kappa       <- exp(x[1])
    nu          <- exp(x[2])
    sigma       <- exp(x[3])
    sigma_noise <- exp(x[4])
    Sigma <- matern.covariance(D, kappa, nu, sigma)
    diag(Sigma)<- diag(Sigma) + sigma_noise^2
    R <- chol(Sigma, pivot=T)
    if(attr(R,"rank" )< length(y))
        return(-Inf)
    v <- solve(t(R),y[attr(R,"pivot")])
    return( -sum(log(diag(R))) - 0.5*t(v)%*%v)
}


#' @param death_prediction - result from uk.prediction
#'                           $date          -
#'                           $date.reported - should only be one date
#'                           $median
#'                           $low (CI)
#'                           $upp (CI)
#'
#' @param   max.days.to.report -  after max.days.to.report data is assume known
####

gp.smooth <- function(death_prediction, max.days.to.report, theta = NULL, CI_width = 0.95){


    if(is.null(theta)){
        theta <- c(0,0,0,0)
    }

    report.date <- unique(death_prediction$date.reported)
    if(length(report.date)>1)
    {
        cat("error only day is allowed in gp.smooth data for input\n")
        return
    }

    model_dt_smooth <- NULL

    index.known <- report.date-death_prediction$date > max.days.to.report

    N <- length(death_prediction$date)
    death_reported_so_far <- data.frame(
                                        Deaths  = death_prediction$median[index.known],
                                        date    = death_prediction$date[index.known])

    time_d  = death_reported_so_far$date - min(death_reported_so_far$date)
    y_d     =  sqrt(death_reported_so_far$Deaths )
    D_d     = as.matrix(dist(time_d))
    res     = optim(theta, function(x){-likelihood_death(x, D_d, y_d)})

    ###
    # approx model with measumerent error
    # kriging
    ###
    y_D <- sqrt(death_prediction$median)
    sigma_y <- (sqrt(death_prediction$upp+1)-sqrt(death_prediction$low+1))/(2*qnorm(0.5+CI_width/2))
    index.obs = is.na(sigma_y)==F

    D <- as.matrix(dist(death_prediction$date-min(death_prediction$date)))
    kappa       <- exp(res$par[1])
    nu          <- exp(res$par[2])
    sigma       <- exp(res$par[3])
    sigma_noise <- exp(res$par[4])
    Sigma <- matern.covariance(D, kappa, nu, sigma)
    Sigma_obs <- Sigma
    diag(Sigma_obs)<- diag(Sigma_obs) + sigma_noise^2
    Sigma_Y <- Sigma_obs
    diag(Sigma_Y)[index.obs] <- diag(Sigma_Y)[index.obs] + sigma_y[index.obs]^2
    mu <- Sigma[,index.obs]%*%solve(Sigma_Y[index.obs,index.obs], y_D[index.obs])
    mu_obs <- Sigma_obs[,index.obs]%*%solve(Sigma_Y[index.obs,index.obs], y_D[index.obs])
    Sigma_cond <- Sigma - Sigma[,index.obs]%*%solve(Sigma_Y[index.obs,index.obs],Sigma[index.obs,])
    Sigma_obs_cond <- Sigma_obs - Sigma_obs[,index.obs]%*%solve(Sigma_Y[index.obs,index.obs],Sigma_obs[index.obs,])

    lw <- ceiling((mu_obs+ qnorm(0.5+CI_width/2)*sqrt(diag(Sigma_obs_cond) + 1e-8))^2)
    uw <- floor(apply(mu_obs- qnorm(0.5+CI_width/2)*sqrt(diag(Sigma_obs_cond) + 1e-8 ),1,function(x){max(0,x)})^2)
    death_prediction$mean <- NA
    death_prediction$median <- round((mu_obs)^2)
    death_prediction$low <- lw
    death_prediction$upp <- uw
    return(death_prediction)
}
