##
# simple stan code for data
#
##
library(rstan)

sim <- 2000
lag <- 0:10
#rstan_options(auto_write = TRUE)

deaths <- readRDS(file.path("data", "processed", "processed_data.rds"))
icu    <- readRDS(file.path("data", "processed", "processed_data_icu.rds"))
icu_pre <-read_excel(path = file.path("data", "FHM", "Folkhalsomyndigheten_Covid19.xlsx"),
           sheet = 3,
           col_types = c("date", "numeric"),
           progress = FALSE)
icu_pre <- data.frame(date = as.Date(icu_pre$Datum_vårdstart),
                      icu  = icu_pre$Antal_intensivvårdade)


Reported = icu$detected
N <- dim(Reported)[1]
incomplete <- 10 # remove the ten latest
icu_est <- apply(Reported, 1, max, na.rm=T)
icu_est <- c(icu_pre$icu[1:51], icu_est[1:(N-incomplete)])
dates <- c(icu_pre$date[1:51] ,as.Date(icu$dates[1:(N-incomplete)]))
data.list <- list(N=length(icu_est),
                  Y = icu_est)



stan_out = stan(file = "src/model2/stan/RW1_negbin.stan",
                data= data.list,
                iter = sim,
                warmup = ceiling(0.5*sim),
                chains =2)
Theta_stan <- rstan::extract(stan_out)$theta
#tau        <- extract(stan_out)$tau

CI <- apply(exp(Theta_stan),2,quantile,c(0.1,0.5,0.9))
plot(dates,icu_est,xlab='days',ylab='icu')

lines(dates,CI[2,],col='blue',lty=2)
lines(dates,CI[1,],col='red',lty=2)
lines(dates,CI[3,],col='red',lty=2)

Cov_analysis <- data.frame(dates = deaths$dates[deaths$dates%in%dates],
                           deaths = apply(deaths$detected[deaths$dates%in%dates,],1,max, na.rm=T))
#

med_ <- apply(Theta_stan,2,quantile,c(0.5))
icu_cov <- c()

for(i in lag){
    icu_lag <- c(rep(0,i), med_)[1:length(med_)]
    icu_cov <- cbind(icu_cov, icu_lag)
}
colnames(icu_cov) <- paste("lag_",lag,sep="")
icu_cov <- icu_cov[dates%in%Cov_analysis$dates,]
lik <- rep(0,length(lag))
for(i in 1:length(lag)){
    X <- icu_cov[,i]
    index <- X>0
    X <- X[index]
    y <- Cov_analysis$deaths[index]
    glm.fit <- glm(y~x,poisson, data= data.frame(y=Cov_analysis,x=X))
    lik[i] <-  logLik(glm.fit)[1]/sum(index)
}
x11()
plot(lag,lik)
i <- which.max(lik)
index <- X>0
X <- X[index]
y <- Cov_analysis$deaths[index]
glm.fit <- glm(y~x,poisson, data= data.frame(y=Cov_analysis,x=X))
lik[i] <-  logLik(glm.fit)[1]/sum(index)
icu    <- readRDS(file.path("data", "processed", "processed_data_icu.rds"))
lag_est <- rep(0, N - incomplete)
for(j in 1:(dim(icu$detected)[1] -incomplete)){
    lag_est[j] <- max(icu$detected[j,1:(j+7)],na.rm=T)
}

icu_cov <- data.frame(date= dates[dates%in%Cov_analysis$dates], icu = X)
saveRDS(icu_cov, file.path("data", "processed", "icu_covariates.rds"))
