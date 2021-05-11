##
# UK analysis
# examning what covariates and how the probabilites works for uk
#
##

source(file.path("src","util","util_swe.r"))
result <- readRDS(file.path("data", "processed", "processed_data_sweden.rds"))


#remove data before  2020-07-01
index <- result$dates_report >= as.Date("2020-07-01")
result$detected             = result$detected[index,index]
result$report               = result$report[index,index]
result$dates                = result$dates[index]
result$dates_report         = result$dates_report[index]
result$dates_not_reported   = result$dates_not_reported[index]
max.days.to.report <- 30



N <- deaths_at_t(result$detected, max.days.to.report)
fixed.dates <- sum(is.na(N)==F)
N.fixed <- N[1:fixed.dates]
Data <- newDeaths(N.fixed,
                  reports = result$detected[1:fixed.dates,1:fixed.dates ],
                  maxusage.day = max.days.to.report)
reports <- result$report[1:fixed.dates,1:fixed.dates]
Data$report.new[ is.na(reports)==F &reports==0 ] = NA
dates_report <- result$dates_report[1:fixed.dates]
X_T <- X.swe(dates_report, reports)
beta <- fit.mu.M(result, max.days.to.report, zero.inflation=T, use.reports = F)
X_T <- X_T[,beta$index ]
mu1 <- 1/(1+exp(-X_T%*%beta$mu1))
M1  <- exp(X_T%*%beta$M1)
alpha1 <- mu1*M1
beta1  <- M1-alpha1



X.zero <- X.swe.zero(dates_report)
pi =1/(1+exp(-X.zero%*%beta$pi))

Alpha1 <- matrix(NA,
                 ncol=fixed.dates,
                 nrow= fixed.dates)
Beta1  <- Pi <- Alpha1

Alpha1[upper.tri(Data$report.new,diag=T)] <- alpha1
Beta1[upper.tri(Data$report.new,diag=T)]  <- beta1
Pi[upper.tri(Data$report.new,diag=T)]     <- pi



lag.plot = 2
MS.plot = F
lag.data.simulated <- zero.BB.dist.by.lag(Data,Alpha1, Beta1, Pi , lags = 1:10)

data.lag <- data.by.lag(Data, lags = 1:10)
days <- c("Sunday")
data.lag$MS         = weekdays(data.lag$date.reported)%in%days
lag.data.simulated$MS  =  weekdays(lag.data.simulated$date.reported)%in%days

plot(data.lag[lag==lag.plot & MS == MS.plot,date.reported], data.lag[lag==lag.plot  & MS == MS.plot, y/n], pch=2,ylim = c(-0.1,1), xlim = c(min(data.lag[,date.reported]),max(data.lag[,date.reported])))
points(lag.data.simulated[lag==lag.plot & MS == MS.plot ,date.reported], lag.data.simulated[lag==lag.plot  & MS == MS.plot,upp/n] ,col='blue')
points(lag.data.simulated[lag==lag.plot  & MS == MS.plot,date.reported], lag.data.simulated[lag==lag.plot  & MS == MS.plot,low/n] ,col='red')


