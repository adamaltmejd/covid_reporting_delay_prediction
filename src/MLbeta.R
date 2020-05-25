
source("src/util.r")
##
#'
#'    lag    - (int)   lag days
#'    theta0 - (p x 1) inital guess of parameters
##
ML_betaBin <- function(result,
                       lag,
                       theta0 = NULL){

    report <- splitlag(result$detected, lag)
    report_lag <- report$Reported_T
    deaths_est_T <- apply(report$Reported_T, 1, max, na.rm=T)

    data_T <- newDeaths(deaths_est_T,
                        report$Reported_T)
    N <- dim(data_T$death.remain)[2]
    X <- setup_data_lag(N, 2, result$dates_report[1:N], 2)
    IndexX <-  row(data_T$report.new) <= (N-lag)
    X <- X[IndexX[upper.tri(data_T$death.remain,diag = T)],]
    index <- upper.tri(data_T$death.remain,diag = T) & row(data_T$report.new) <= (N-lag)
    y = data_T$report.new[index ]
    n = data_T$death.remain[index ]
    index = is.na(y)==F & is.na(n)==F
    X <- as.matrix(X[index,])
    y <- y[index]
    n <- n[index]
    p <- dim(X)[2]
    if(is.null(theta0))
        theta0 <- rep(0,2*dim(X)[2])
    res <- optim(theta0, method='CG',function(x){-log_bb(x, y, n, X)})
    return(list(alpha_X = res$par[1:p], beta_X = res$par[(p+1):(2*p)]))
}

log_bb<- function(theta, y, n,X){
    p <- length(theta)/2
    alpha <- exp(X%*%theta[1:p])
    beta  <- exp(X%*%theta[(p+1):(2*p)])
    sum(dBB(x =y,
            size = n,
            alpha = alpha,
            beta = beta,
            log.p=T))
}
##
#' running the benchmark for day j
#'
#' j - day j
#' result    - (output) of data
#' true.day  - (int) days not estimated
#'  lag    - (int)   lag days
#' MCMC_sim     - (int) how many simulations
#' burnin_p     - ([0,1]) how many burnins to run in percantage of MCMC_sim
#' prior        - (1 ) [1] 0: rw1, 1: rw2,
##
benchmark_BetaGP_lag_j <-function(j,
                                  result,
                                  lag,
                                  MCMC_sim,
                                  burnin_p,
                                  deaths_sim,
                                  prior = c(0)){
    require(Matrix)
    N_T <- length(result$dates)
    deaths_est_T <- apply(result$detected, 1, max, na.rm=T)
    result_j <- result
    result_j$detected     <- result_j$detected[1:j,1:j]
    result_j$dates        <- result_j$dates[1:j]
    result_j$dates_report <-  result_j$dates[1:j]
    report_j <- splitlag(result_j$detected, lag)
    param<- ML_betaBin(result_j,
                       lag)
    N_j <- j
    X_T <- setup_data_lag(N_T, 2, result$dates_report, 2)
    X_j <- setup_data_lag(N_j, 2, result_j$dates_report, 2)


    deaths_est <- apply(report_j$Reported_T, 1, max, na.rm=T)
    ##
    # seting up MCMC GP
    ##
    Start_theta <- log(deaths_est+1)
    tau <- 0.1
    MH_obj_GP <- MH_setup()
    MH_obj_GP$sigma <- 0.1

    MH_obj_GP$theta <- Start_theta
    if(prior == 0){
        n_theta <- length(Start_theta)
        L <- toeplitz(c(-1,1, rep(0,n_theta-2)))
        L[lower.tri(L,diag=F)] <- 0
        L <- L[-n_theta,]
        L<-as(L, "sparseMatrix")
    }else{
        n_theta <- length(Start_theta)
        L <- toeplitz(c(-1,2,-1, rep(0,n_theta-3)))
        L[lower.tri(L,diag=F)] <- 0
        L <- L[-c(n_theta,n_theta-1),]
        L<-as(L, "sparseMatrix")
    }
    MH_obj_GP$Lik <- function(x, N, L) {
        prior_list <- prior_GP(x, L)
        lik_list   <- lik_poisson(N, x)
        return(list(loglik = lik_list$loglik + prior_list$loglik,
                    grad = lik_list$grad  + prior_list$grad,
                    Hessian = lik_list$Hessian + prior_list$Hessian))
    }


    Alpha <- matrix(NA, ncol=N_j, nrow=N_j)
    Beta <- matrix(NA, ncol=N_j, nrow=N_j)
    P <- matrix(NA, ncol=N_j, nrow=N_j)
    theta_GP <-  matrix(NA, nrow=MCMC_sim, ncol = length(MH_obj_GP$theta))
    GPprior <- rep(0, MCMC_sim)
    Death_est <- matrix(NA, nrow=MCMC_sim, ncol=N_j)
    alpha.MCMC <- rep(4, N_j)
    pred_set <-array(NA, dim = c(j,N_T-j,MCMC_sim))
    p <- dim(X)[2]

    burnin <- ceiling(burnin_p*MCMC_sim)


    data_ <-newDeaths(deaths_est,
                      report_j$Reported_T,
                      Inf)

    for(i in 1:(MCMC_sim+burnin-1)){

        if(i%%100==0){
            cat('*')
        }
        if(i%%1000==0){
            cat('+')
        }
        if(i%%10000==0){
            cat(' ',i/10000,' ')
        }

        Alpha[upper.tri(data_$report.new,diag=T)] <- exp(X_j%*%param$alpha_X)
        Beta[upper.tri(data_$report.new,diag=T)]  <-  exp(X_j%*%param$beta_X)

        prior_N <- function(N, i){ N *MH_obj_GP$theta[i] - lgamma(N+1)}
        res <- sample.deathsBB(deaths_sim,
                               deaths_est,
                               Alpha,
                               Beta,
                               report_j$Reported_T,
                               alpha.MCMC,
                               true.day = N_j - lag,
                               use.prior=T,
                               prior=prior_N)
        deaths_est <- res$deaths
        data_ <-newDeaths(deaths_est,
                          report_j$Reported_T,
                          Inf)

        ##
        # build tau sampler
        ##

        L_theta          <- L%*%MH_obj_GP$theta
        tau   <- rgamma(1, shape=length(L_theta)/2 + 0.001, scale = sum(L_theta^2)/2 + 0.001)
        L_in <- sqrt(tau) * L

        MH_obj_GP <- MALAiter(MH_obj_GP, TRUE,
                              N = deaths_est,
                              L = L_in)

        if(i < burnin){
            alpha.MCMC[res$acc/deaths_sim > 0.3] <- alpha.MCMC[res$acc/deaths_sim > 0.3] +1
            alpha.MCMC[res$acc/deaths_sim < 0.3] <- alpha.MCMC[res$acc/deaths_sim < 0.3] -1
            alpha.MCMC[alpha.MCMC<1] <- 1
        }else{
            Death_est[i-burnin + 1,] <-res$deaths
            theta_GP[i-burnin + 1,] <-  as.vector(MH_obj_GP$theta)
            GPprior[i-burnin + 1]  <- tau
        }
    }
    res_save <- list(Death_est = Death_est,
                     theta_GP = theta_GP,
                     lag       = lag,
                     GPprior      = GPprior,
                     date  = result$dates_report[j],
                     j = j)
    return(res_save)

}
