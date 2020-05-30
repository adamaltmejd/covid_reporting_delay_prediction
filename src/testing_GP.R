###
# GP test
#
#
###
#kappa  =0.3
source(file.path("src", "GPutil.R"))

result <- readRDS(file.path("data", "processed", "processed_data.rds"))
nu <- 1.8
Reported = result$detected
N <- dim(Reported)[1]

deaths_est <- apply(Reported, 1, max, na.rm=T)
rem <- 20
death_est <- deaths_est[1:(N-rem)]
loc <- 1:(N-rem)

h <- as.matrix(dist(loc, method = "euclidean"))
One <- rep(1, length(loc))
OneOneT <- One%*%t(One)
#kappa, nu, sigma, sigma_e
theta <- c(0,0,0,0)
f <-  function(x){-lik.spatial(log(death_est),matern.covariance(h,
                                                                exp(x[1]),
                                                                nu,
                                                                exp(x[2])) +
                                   exp(x[3])*diag(length(loc)) + 10^3*OneOneT)}
res <- optim(theta,f)
theta <- exp(res$par)
kappa <- theta[1]
sigma <- theta[2]
sigma_e <- sqrt(theta[3])
Cov <- matern.covariance(h,
                         kappa,
                         nu,
                         sigma) +
    sigma_e^2*diag(length(loc)) + 10^3*OneOneT
n <- length(loc)
pred_loc <- 20
loc <- c(loc,loc[n]+1:(pred_loc))

h_pred <- as.matrix(dist(loc, method = "euclidean"))


Sigma_11 <- matern.covariance(h_pred,
                              kappa,
                              nu,
                              sigma) +
    + 10^3*rep(1, length(loc))%*%t(rep(1, length(loc)))
post_mean <- Sigma_11[,1:n]%*%solve(Cov,log(death_est))
post_var  <- Sigma_11 - Sigma_11[,1:n]%*%solve(Cov,t(Sigma_11[,1:n]))
plot((death_est),xlim=c(0,n+pred_loc),main=format(result$dates_report[N-rem]))
lines(exp(post_mean))
lines(exp(post_mean+qnorm(0.05)*sqrt(diag(post_var))),col='red')
lines(exp(post_mean+qnorm(0.95)*sqrt(diag(post_var))),col='red')
print(kappa)
