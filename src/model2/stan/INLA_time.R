# simple stan code for data
#
##
library(INLA)

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
icu_est <- rep(0,N-incomplete)
for(j in 1:(N -incomplete)){
    icu_est[j] <- max(icu$detected[j,1:(j+8)],na.rm=T)
}


icu_est <- c(icu_pre$icu[1:51], icu_est)
dates <- c(icu_pre$date[1:51] ,as.Date(icu$dates[1:(N-incomplete)]))
theta <- rep(0,length(icu_est))
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
}

plot(dates,icu_est,xlab='days',ylab='icu')
lines(dates,exp(theta))
icu_cov <- c()

for(i in lag){
    icu_lag <- c(rep(0,i), theta)[1:length(theta)]
    icu_cov <- cbind(icu_cov, icu_lag)
}
colnames(icu_cov) <- paste("lag_",lag,sep="")
icu_cov <- icu_cov[dates%in%Cov_analysis$dates,]
lik <- rep(0,length(lag))
for(i in 1:length(lag)){
    X <- icu_cov[,i]
    y <- Cov_analysis$deaths
    glm.fit <- glm(y~x,poisson, data= data.frame(y=y,x=X))
    lik[i] <-  logLik(glm.fit)[1]/sum(index)
}
x11()
w <- exp(lik-max(lik))/sum(exp(lik-max(lik)))
plot(lag,w)

i <- 8
X <- icu_cov[,i]
index <- X>0
X <- X[index]
y <- Cov_analysis$deaths[index]
glm.fit <- glm(y~x,poisson, data= data.frame(y=y,x=X))


icu_cov <- data.frame(date= dates[dates%in%Cov_analysis$dates], icu = X)
saveRDS(icu_cov, file.path("data", "processed", "icu_covariates.rds"))

