library(rstan)


source("src/util.r")
sim <- 1000
unique.days = 5
maxusage.day=20
result <- readRDS(file.path("data", "processed", "processed_data.rds"))
Reported_T = result$detected
N_T <- dim(Reported_T)[1]
fixed = N_T

holidays.Sweden <- as.Date(c("2020-04-10","2020-04-13","2020-05-01","2020-05-21"))
holidays <- weekdays(result$dates_report)%in%c("Sunday","Saturday") | c(result$dates_report)%in%c(holidays.Sweden)
Reported_T = result$detected[1:fixed,1:fixed]
Reported_O = result$detected[1:fixed,1:fixed]
N_T <- dim(Reported_T)[2]
lag <- 3
for(i in 1:N_T){
    k <- which.max(cumsum(holidays[i:N_T]==F)==(lag-1))
    Reported_O[i, i:min(i+k, N_T)] <- NA
    if(k+i<N_T){

        Reported_O[i,  (i+k+1):N_T] <- Reported_O[i,  (i+k+1):N_T] - Reported_T[i,  (i+k)]
        Reported_T[i, (i+k+1):N_T]<- NA

    }

}

deaths_est_T <- apply(Reported_T, 1, max, na.rm=T)

data_T <- newDeaths(deaths_est_T,
                    Reported_T,
                    maxusage.day =maxusage.day)


deaths_est_O <- apply(Reported_O, 1, max, na.rm=T)
data_O <- newDeaths(deaths_est_O,
                    Reported_O,
                    maxusage.day =maxusage.day)
X_O <- setup_data(fixed, maxusage.day, result$dates_report[1:fixed], 2)
index <- upper.tri(data_O$death.remain,diag = T)
y_O = data_O$report.new[index]
n_O = data_O$death.remain[index]
index = is.na(y)==F & is.na(n)==F
X_O <- as.matrix(X_O[index,])
X_O <- X_O[,colSums(X_O)>0]
y_O <- y_O[index]
n_O <- n_O[index]
fit_O <- glm( cbind(y_O, n_O) ~ X_O,
                 family = "binomial")

X_T <- setup_data(fixed, maxusage.day, result$dates_report[1:fixed], 2)
index <- upper.tri(data_T$death.remain,diag = T)
y_T = data_T$report.new[index]
n_T = data_T$death.remain[index]
index = is.na(y)==F & is.na(n)==F
X_T <- as.matrix(X_T[index,])
X_T <- X_T[,colSums(X_T)>0]
y_T <- y_T[index]
n_T <- n_T[index]
fit_T <- glm( cbind(y_T, n_T) ~ X_T,
                 family = "binomial")
print(summary(fit_O))
print(summary(fit_T))
par(mfrow=c(1,2))
plot(deaths_est_T)
plot(deaths_est_O)
