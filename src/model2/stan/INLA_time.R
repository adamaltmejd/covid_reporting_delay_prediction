# simple stan code for data
#
##
library(mgcv)
lag <- 0:14
#rstan_options(auto_write = TRUE)

deaths <- readRDS(file.path("data", "processed", "processed_data.rds"))
icu    <- readRDS(file.path("data", "processed", "processed_data_icu.rds"))
icu_pre <-read_excel(path = file.path("data", "Folkhalsomyndigheten_Covid19.xlsx"),
                     sheet = 3,
                     col_types = c("date", "numeric"),
                     progress = FALSE)
icu_pre <- data.frame(date = as.Date(icu_pre$Datum_vårdstart),
                      icu  = icu_pre$Antal_intensivvårdade)



Reported = icu$detected
N <- dim(Reported)[1]
incomplete <- 10 # remove the ten latest
incomplete_deaths <- 30
icu_est <- rep(0,N-incomplete)
for(j in 1:(N -incomplete)){
    icu_est[j] <- max(icu$detected[j,1:(j+8)],na.rm=T)
}


icu_est <- c(icu_pre$icu[1:51], icu_est)
dates <- c(icu_pre$date[1:51] ,as.Date(icu$dates[1:(N-incomplete)]))
theta <- theta2 <- rep(0,length(icu_est))
i <- 51
data.list <- list(date=1:i,
                  Y = icu_est[1:i])
model.fit <- gam(Y~s(date),data=data.listfamily  =poisson)
theta[1:i] <- log(model.fit$fitted.values)[1:i]
for(i in 52:length(icu_est)){
    print(i)
    data.list <- list(date=1:i,
                  Y = icu_est[1:i])
    model.fit <- gam(Y~s(date),data=data.list,family= poisson)
    theta[i] <- log(model.fit$fitted.values)[i]
}

plot(dates,icu_est,xlab='days',ylab='icu')
lines(dates,exp(theta),col='blue')
icu_cov  <- c()
for(i in lag){
    icu_lag <- c(rep(0,i), theta)[1:length(theta)]
    icu_cov <- cbind(icu_cov, icu_lag)
}
colnames(icu_cov) <- paste("lag_",lag,sep="")
N_d <- length(deaths$dates)
Cov_analysis <- data.frame(deaths = apply(deaths$detected[1:(N_d-incomplete_deaths),],1,max,na.rm=T),
                           dates = deaths$dates[1:(N_d-incomplete_deaths)])
icu_cov <- icu_cov[dates%in%Cov_analysis$dates,]
lik  <- rep(0,length(lag))
y <- Cov_analysis$deaths

for(i in 1:length(lag)){
    X <- icu_cov[,i]
    glm.fit <- glm(y~x,poisson, data= data.frame(y=y,x=X))
    lik[i] <-  logLik(glm.fit)[1]/length(y)

}
x11()
lik_joint <- c(lik)
w1 <- exp(lik-max(lik_joint))/sum(exp(lik-max(lik_joint)))
plot(lag,w1,ylim=c(0,0.2),main='which lag')
x11()
plot(Cov_analysis$dates,Cov_analysis$deaths)

i <- 8
X <- icu_cov[,i]
glm.fit <- glm(y~x,poisson, data= data.frame(y=y,x=X))
lines(Cov_analysis$dates,exp(cbind(1,X)%*%glm.fit$coefficients),col='red',lty=3)

icu_cov <- data.frame(date= dates[dates%in%Cov_analysis$dates], icu = X)
saveRDS(icu_cov, file.path("data", "processed", "icu_covariates.rds"))

