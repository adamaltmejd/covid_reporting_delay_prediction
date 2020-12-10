
library(INLA)
data <- readRDS(file.path("data", "processed", "processed_data.rds"))
N <- dim(data$detected)[1]
icu_cov  <- readRDS(file.path("data", "processed", "icu_covariates.rds"))
deaths <- apply(data$detected,1,max,na.rm=T)[1:(N-30)]
dates <- data$dates[1:(N-30)]
data.list <- list(date=1:length(dates),
                  icu  =  icu_cov$icu[icu_cov$date%in%dates],
                  Y = deaths)

model <- inla(
    Y ~ icu +
        f(date, model = "iid"),
    data = data.list,
    family = "poisson"
)

sigma <- sqrt(1/model$summary.hyperpar$"0.5quant")
mode <-  model$summary.fixed$"0.5quant"
#mode <- c(0.9139627, 0.9819359)
print(mode)
sim = 1000


plot(dates, deaths, xlim = c(min(icu_cov$date),max(icu_cov$date)))
lines(icu_cov$date,exp(cbind(1,icu_cov$icu)%*%mode),col='red')
Quan_Y <- matrix(0, nrow=length(icu_cov$date), ncol=2)
for(i in 1:length(icu_cov$date)){
    Z = sigma*rnorm(sim)
    Ys <- rpois(sim, exp(mode[1]+icu_cov$icu[i]*mode[2]+Z))
    Quan_Y[i,] = quantile(Ys,c(0.025,0.975))
}
lines(icu_cov$date,Quan_Y[,1],col='blue')
lines(icu_cov$date,Quan_Y[,2],col='blue')
cat('lower = ',mean(Quan_Y[,1]>deaths),'\n')
cat('upper = ',mean(Quan_Y[,2]<deaths),'\n')
