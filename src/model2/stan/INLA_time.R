# simple stan code for data
#
##
library(INLA)

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

model <- inla(
    Y ~
        f(date, model = "rw2"),
    data = data.list,
    family = "poisson"
)
theta[1:i] <- model$summary.fixed$"0.5quant"+model$summary.random$date$"0.5quant"[1:i]
theta2[1:i] <- model$summary.fixed$"0.5quant"+model$summary.random$date$"0.5quant"[1:i]
for(i in 52:length(icu_est)){
    print(i)
    data.list <- list(date=1:i,
                  Y = icu_est[1:i])

    model <- inla(
        Y ~
            f(date, model = "rw2"),
        data = data.list,
        family = "poisson"
    )
    theta[i] <- model$summary.fixed$"0.5quant"+model$summary.random$date$"0.5quant"[i]
    theta2[i] <- model$summary.fixed$"0.5quant"+model$summary.random$date$"0.5quant"[i-1]
}

plot(dates,icu_est,xlab='days',ylab='icu')
lines(dates,exp(theta),col='blue')
lines(dates,exp(theta2),col='red')
icu_cov  <- c()
icu_cov2 <- c()
for(i in lag){
    icu_lag <- c(rep(0,i), theta)[1:length(theta)]
    icu_cov <- cbind(icu_cov, icu_lag)
    icu_lag2 <- c(rep(0,i), theta2)[1:length(theta)]
    icu_cov2 <- cbind(icu_cov2, icu_lag2)
}
colnames(icu_cov) <- paste("lag_",lag,sep="")
colnames(icu_cov2) <- paste("lag_",lag,sep="")
N_d <- length(deaths$dates)
Cov_analysis <- data.frame(deaths = apply(deaths$detected[1:(N_d-incomplete_deaths),],1,max,na.rm=T),
                           dates = deaths$dates[1:(N_d-incomplete_deaths)])
icu_cov <- icu_cov[dates%in%Cov_analysis$dates,]
icu_cov2 <- icu_cov2[dates%in%Cov_analysis$dates,]
lik  <- rep(0,length(lag))
lik2 <- rep(0,length(lag))
y <- Cov_analysis$deaths

for(i in 1:length(lag)){
    X <- icu_cov[,i]
    glm.fit <- glm(y~x,poisson, data= data.frame(y=y,x=X))
    lik[i] <-  logLik(glm.fit)[1]/sum(index)

    X2 <- icu_cov2[,i]
    glm.fit <- glm(y~x,poisson, data= data.frame(y=y,x=X2))
    lik2[i] <-  logLik(glm.fit)[1]/sum(index)

}
x11()
lik_joint <- c(lik,lik2)
w1 <- exp(lik-max(lik_joint))/sum(exp(lik-max(lik_joint)))
w2 <- exp(lik2-max(lik_joint))/sum(exp(lik2-max(lik_joint)))
plot(lag,w1,ylim=c(0,0.2),main='which lag')
points(lag, w2,col='red',pch=4)
x11()
plot(Cov_analysis$dates,Cov_analysis$deaths)

i <- 8
X <- icu_cov[,i]
X2 <- icu_cov2[,i]
glm.fit <- glm(y~x,poisson, data= data.frame(y=y,x=X))
glm.fit2 <- glm(y~x,poisson, data= data.frame(y=y,x=X2))
lines(Cov_analysis$dates,exp(cbind(1,X)%*%glm.fit$coefficients),col='red',lty=3)
lines(Cov_analysis$dates,exp(cbind(1,X2)%*%glm.fit$coefficients),col='blue',lty=3)

icu_cov <- data.frame(date= dates[dates%in%Cov_analysis$dates], icu = X)
saveRDS(icu_cov, file.path("data", "processed", "icu_covariates.rds"))

