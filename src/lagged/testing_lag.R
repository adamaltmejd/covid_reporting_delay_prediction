library(rstan)


source("src/util.r")
source("src/MLbeta.R")
sim <- 1000
unique.days = 5
maxusage.day=20
result <- readRDS(file.path("data", "processed", "processed_data.rds"))
Reported_T = result$detected
N_T <- dim(Reported_T)[1]
fixed = N_T
lag <- 3

report <- splitlag(result$detected[1:fixed,1:fixed],as.Date(result$dates_report[1:fixed]), lag)
deaths_O_est <- apply(report$Reported_O, 1, max, na.rm=T)

data_O <-newDeaths(deaths_O_est,
                  report$Reported_O,
                  maxusage.day)
X_O <- setup_data_postlag(fixed, 2, result$dates_report[1:fixed], 1)
X_O <- X_O[,colSums(X_O)>0]
index <- upper.tri(data_O$death.remain,diag = T)
y_O = data_O$report.new[index]
n_O = data_O$death.remain[index]
index = is.na(y_O)==F & is.na(n_O)==F
X_O <- as.matrix(X_O[index,])

y_O <- y_O[index]
n_O <- n_O[index]
fit_O <- glm( cbind(y_O, n_O) ~ X_O-1,
                 family = "binomial")

res1 <- optim(rep(0, dim(X_O)[2]*2),
             function(x){-loglikProbBB(x, data_O$death.remain, data_O$report.new, X_O)$loglik},
             function(x){-loglikProbBB(x, data_O$death.remain, data_O$report.new, X_O)$grad},
             method = 'CG')
res <- optim(res1$par,
             function(x){-loglikProbBB(x, data_O$death.remain, data_O$report.new, X_O)$loglik})



###
# XT
#
###
X_T <- setup_data_postlag(fixed, maxusage.day, result$dates_report[1:fixed], 2)
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
