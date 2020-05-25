library(rstan)


source("src/util.r")
sim <- 1000
unique.days = 5
maxusage.day=20
result <- readRDS(file.path("data", "processed", "processed_data.rds"))
Reported_T = result$detected
N_T <- dim(Reported_T)[1]
fixed = N_T
Reported_T <- result$detected[1:fixed,1:fixed]
deaths_est_T <- apply(Reported_T, 1, max, na.rm=T)
deaths_est_T <- deaths_est_T[1:fixed]
data_T <- newDeaths(deaths_est_T,
                    Reported_T[1:fixed,1:fixed],
                    maxusage.day =maxusage.day)
death.remain = data_T$death.remain
report.new   = data_T$report.new
X_T <- setup_data(fixed, maxusage.day, result$dates_report[1:fixed], unique.days)

index <- upper.tri(death.remain,diag = T)
y = report.new[index]
n = death.remain[index]
index = is.na(y)==F
X <- as.matrix(X_T[index,])
y <- y[index]
n <- n[index]
data.list<- list(Y=y,
                 n = n,
                 N = length(n),
                 P = dim(X)[2],
                 X = X,
                 sigma_beta = 3,
                 sigma_alpha = 3)
stan_out = stan(file = "src/stan/betabin.stan",
                data= data.list,
                iter = sim,
                warmup = ceiling(0.5*sim),
                chains =2,
                cores=1)

alpha <- t(exp(X%*%t(extract(stan_out)$beta_alpha)))
beta <- t(exp(X%*%t(extract(stan_out)$beta_beta)))
P_obs <-alpha/(alpha+beta)
i <- 1
obs <- y[X[,i]==1]/n[X[,i]==1]
print(colMeans((P_obs[,X[,i]==1])))
index<-colMeans((P_obs[,X[,i]==1]))==colMeans((P_obs[,X[,i]==1]))[1]
plot(colMeans((P_obs[,X[,i]==1])),obs)
print(mean(obs[index], na.rm=T))

alpha_2 <- t(exp(X%*%t(res_save$Thetas[,1:10])))
beta_2 <- t(exp(X%*%t(res_save$Thetas[,11:20])))
P_obs2 <-alpha_2/beta_2
plot(colMeans((P_obs[,X[,i]==1])),colMeans((P_obs2[,X[,i]==1])))

k <- 599
j <- 389
alpha <- exp(X%*%extract(stan_out)$beta_alpha[k,])
beta <-  exp(X%*%extract(stan_out)$beta_beta[k,])

alpha_mod <- exp(X%*%res_save$Thetas[k,1:10])
beta_mod <- exp(X%*%res_save$Thetas[k,11:20])
y_1 <- y[j]
n_1 <- n[j]
alpha_1 <- alpha[j]
beta_1 <- beta[j]
alpha_1s <- alpha_mod[j]
beta_1s <- beta_mod[j]
par(mfrow=c(1,2))
plot(1:100,dBB(1:100,n_1,alpha_1,beta_1,log.p = F),type='l',
     main=paste('Y=',y_1))
abline(v=y_1,col='green')
lines(1:100,dBB(1:100,n_1,alpha_1s,beta_1s,log.p = F),col='red')

plot(1:150,dBB(y_1,1:150,alpha_1,beta_1,log.p = F),type='l',
     main=paste('n=',y_1))
abline(v=n_1,col='green')
lines(1:150,dBB(y_1,1:150,alpha_1s,beta_1s,log.p = F),col='red')
print(which(X[,2]==1))

load(files[i+2])
lag <- 8
i <- 34
k <- 10
X_Base <- setup_data(N_T, maxusage.day, result$dates_report, unique.days)
Alpha_T <- matrix(NA, N_T,N_T)
Beta_T  <- matrix(NA, N_T,N_T)
Alpha_T2 <- matrix(NA, N_T,N_T)
Beta_T2  <- matrix(NA, N_T,N_T)

Alpha_T[upper.tri(Alpha_T,diag=T)] <- exp(X_Base%*%extract(stan_out)$beta_alpha[k,])
Beta_T[upper.tri(Beta_T,diag=T)]   <- exp(X_Base%*%extract(stan_out)$beta_beta[k,])
Alpha_T2[upper.tri(Alpha_T,diag=T)] <- exp(X_Base%*%res_save$Thetas[k,1:10])
Beta_T2[upper.tri(Beta_T,diag=T)]   <- exp(X_Base%*%res_save$Thetas[k,11:20])
Reported <- result$detected[1:N_T,1:N_T]
Reported_i <- Reported[i,i:(i+lag)]
alpha_i    <- Alpha_T[i,i:(i+lag)]
beta_i     <- Beta_T[i,i:(i+lag)]
alpha_i2    <- Alpha_T2[i,i:(i+lag)]
beta_i2     <- Beta_T2[i,i:(i+lag)]
deaths <-  Reported_i[lag+1]:150
log_d  <- rep(0,length(deaths))
log_d2  <- rep(0,length(deaths))
for(j in 1:length(deaths)){
    log_d[j] <- loglikDeathsGivenProbBB(deaths[j],alpha_i , beta_i, Reported_i)
    log_d2[j] <- loglikDeathsGivenProbBB(deaths[j],alpha_i2 , beta_i2, Reported_i)
}
par(mfrow=c(1,1))
plot(deaths,exp(log_d)/sum(exp(log_d)),type='l',
     main=paste(format(result$dates_report[i]),' Truth=',max(Reported[i,],na.rm=T)))
lines(deaths,exp(log_d2)/sum(exp(log_d2)),col='red')
abline(v=max(Reported[i,],na.rm=T),col='green')
