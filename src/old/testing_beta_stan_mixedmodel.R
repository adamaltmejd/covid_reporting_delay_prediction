library(rstan)
rstan_options(auto_write = TRUE)

source("src/util.r")
sim <- 2000
unique.days = 5
maxusage.day=20

result <- readRDS(file.path("data", "processed", "processed_data.rds"))
Reported_T = result$detected
N_T <- dim(Reported_T)[1]
fixed = 45
Reported_T <- result$detected[1:fixed,1:fixed]
deaths_est_T <- apply(result$detected[1:fixed,], 1, max, na.rm=T)

X_days <- setup_data_mixture_cov(fixed, 2, result$dates_report[1:fixed])
X_trend <- setup_all_days(fixed)/7

X_mixed <- setup_data_mixed_effect(fixed, 2, result$dates_report[1:fixed])
X_mix_obs <- X_days
X_mu <- setup_data(fixed, maxusage.day,  result$dates_report[1:fixed], 1)
X_M  <- X_mu[,1:2]
X_longlag <- rowSums(X_mixed)>0
X_mu <- cbind(X_mu, X_trend* X_longlag, X_trend*(X_longlag==F))
X_M  <- cbind(X_M[,1],X_M[,2]*(X_longlag==0), X_longlag)


data_T <- newDeaths(deaths_est_T,
                    Reported_T,
                    maxusage.day =maxusage.day)
death.remain = data_T$death.remain
report.new   = data_T$report.new


index <- upper.tri(death.remain,diag = T)
y = report.new[index]
n = death.remain[index]
X_mix_obs = X_mix_obs[index]
index = is.na(y)==F & n>0
X_mu <- as.matrix(X_mu[index,])
X_M <- as.matrix(X_M[index,])
X_mix_obs <- X_mix_obs[index]
X_mix_obs_old <- X_mix_obs
X_mix_obs <- as.integer(factor(X_mix_obs))-1

y <- y[index]
n <- n[index]




data.list<- list(Y=y,
                 ndays = max(X_mix_obs),
                 day = X_mix_obs,
                 n = n,
                 N = length(n),
                 P_mu = dim(X_mu)[2],
                 P_M = dim(X_M)[2],
                 X_mu = X_mu,
                 X_M  = X_M,
                 sigma_mu = 3,
                 sigma_M = 3)
stan_out = stan(file = "src/stan/betabin_mixedmodel.stan",
                data= data.list,
                iter = sim,
                warmup = ceiling(0.5*sim),
                chains =2,
                cores=1)
stan_simple = stan(file = "src/stan/betabin_2.stan",
                data= data.list,
                iter = sim,
                warmup = ceiling(0.5*sim),
                chains =2,
                cores=1)

k <- 1
beta_mu <- rstan::extract(stan_out,'beta_mu')$beta_mu[k,]
beta_M  <- rstan::extract(stan_out,'beta_M')$beta_M[k,]
mu_day  <- rstan::extract(stan_out,'mu_day')$mu_day[k,]
M_day  <- rstan::extract(stan_out,'mu_day')$mu_day[k,]
mu = 1/(1+exp(-X_mu%*%beta_mu))
M = exp(X_M%*%beta_M);
mu_ = 1/(1+exp(-X_mu%*%beta_mu - mu_day[21]))
M_ = exp(X_M%*%beta_M + M_day[21]);
#P_full[unique(X_mix_obs_old[X_mix_obs_old>0])] <- summary(stan_out,'mu_day')$summary[,1]
#plot(result$dates_report[1:fixed],P_full)


mu_day <- summary(stan_out,'mu_day')$summary[,1]
names(mu_day) <- result$dates_report[unique(X_mix_obs_old[X_mix_obs_old>0])]

#M_day <- summary(stan_out,'M_day')$summary[,1]
#names(M_day) <- result$dates_report[unique(X_mix_obs_old[X_mix_obs_old>0])]
#print(cbind(mu_day,M_day))
beta_mu <- summary(stan_out,'beta_mu')$summary[,1]
names(beta_mu)<- colnames(X_mu)
print(beta_mu)
beta_M <- summary(stan_out,'beta_M')$summary[,1]
names(beta_M) <- colnames(X_M)
print(beta_M)

beta_mu_ <- summary(stan_simple,'beta_mu')$summary[,1]
names(beta_mu_)<- colnames(X_mu)
print(cbind(beta_mu_,beta_mu))
beta_M_ <- summary(stan_simple,'beta_M')$summary[,1]
names(beta_M_) <- colnames(X_M)
print(cbind(beta_M_,beta_M))
