library(rstan)
rstan_options(auto_write = TRUE)

source("src/util.r")
sim <- 2000
unique.days = 5
maxusage.day=20

result <- readRDS(file.path("data", "processed", "processed_data.rds"))
Reported_T = result$detected
N_T <- dim(Reported_T)[1]
fixed = N_T-15
Reported_T <- result$detected[1:fixed,1:fixed]
deaths_est_T <- apply(result$detected[1:fixed,], 1, max, na.rm=T)

X_days <- setup_data_mixture_cov(fixed, 2, result$dates_report[1:fixed])
X_mix_obs <- X_days
X_mu <- setup_data(fixed, maxusage.day,  result$dates_report[1:fixed], 1)
X_M  <- X_mu[,1:2]
X_longlag <- X_days[upper.tri(death.remain,diag = T)]>0
X_mu <- cbind(X_mu, X_longlag)
X_M  <- cbind(X_M, X_longlag)
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
stan_out = stan(file = "src/stan/betabin_mix.stan",
                data= data.list,
                iter = sim,
                warmup = ceiling(0.5*sim),
                chains =2,
                cores=1)
si <- length(extract(stan_out,'mu_add')$mu_add)
P_T <- matrix(0,
              nrow = data.list$ndays,
              ncol = 2
)
for(k in 1:si){
    theta <- extract(stan_out,'theta')$theta[k,]
    beta_mu <- extract(stan_out,'beta_mu')$beta_mu[k,]
    beta_M  <- extract(stan_out,'beta_M')$beta_M[k,]
    mu_add  <- extract(stan_out,'mu_add')$mu_add[k]
    M_sub  <- extract(stan_out,'M_sub')$M_sub[k]
    P_days <- matrix(0,
                     nrow = data.list$ndays,
                     ncol = 2
                     )
    for(i in 1:data.list$ndays)
        P_days[i,] <- log(theta)

    day <- data.list$day
    Y   <- data.list$Y
    n   <- data.list$n
    mu = 1/(1+exp(-X_mu%*%beta_mu))
    M = exp(X_M%*%beta_M);
    alpha = M * mu;
    mu0 = 1/(1+exp(-X_mu%*%beta_mu - mu_add))
    M0 = exp(X_M%*%beta_M - M_sub);
    alpha_0 = M0 * mu0;
    beta  = (1-mu) * M;
    beta0  = (1-mu0) * M0;
    for(i in 1:data.list$N) {
        if(data.list$day[i]>0){
            P_days[day[i],1] = P_days[day[i],1] + dBB(Y[i],n[i],alpha_0[i],beta0[i]);
            P_days[day[i],2] = P_days[day[i],2] + dBB(Y[i],n[i],alpha[i],beta[i]);
        }
    }
    P_days <- exp(P_days-apply(P_days,1,max))
    P_days <- P_days/rowSums(P_days)
    P_T    <- P_T+ P_days
}
P_T <- P_T/si
P_full <- rep(0,fixed)
P_full[unique(X_mix_obs_old[X_mix_obs_old>0])] <- P_T[,1]
plot(result$dates_report[1:fixed],P_full)


