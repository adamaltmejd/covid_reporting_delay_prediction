
##
#' running the benchmark for day j
#'
#' j - day j
#' result -
#' true.day  - (int) days not estimated
#' maxusage.day - (int) how many reports days to use for each death
#' unique.days  - (int) how many unique reports days to have as parameters
#' MCMC_sim     - (int) how many simulations
#' burnin_p     - ([0,1]) how many burnins to run in percantage of MCMC_sim
#' prior        - (3 x 1) kappa , nu, sigma, sigma_mu
##
benchmark_BetaGP_j <- function(j,
                               result,
                               true.day,
                               maxusage.day,
                               unique.days,
                               MCMC_sim,
                               burnin_p,
                               deaths_sim,
                               theta = c(0.5,2,0.09,100)){
    require(Matrix)
    N_T <- dim(result$detected)[1]
    deaths_est_T <- apply(result$detected, 1, max, na.rm=T)
    Reported <- result$detected[1:j,1:j]
    dates_report <- result$dates_report[1:j]
    N <- dim(Reported)[1]
    deaths_est <- apply(Reported, 1, max, na.rm=T)
    deaths_est[1:true.day] = deaths_est_T[1:(true.day)] #days with known deaths


    data_ <- newDeaths(deaths_est,
                       Reported,
                       maxusage.day =maxusage.day)


    ##
    # building covariate matrix
    ##
    # building covariate matrix
    ##
    X_mu <- setup_data(N, maxusage.day,  result$dates_report[1:N], 1)
    X_M  <- X_mu[,1:2]
    X_mixed <- setup_data_mixed_effect(N, 2, result$dates_report[1:N])
    X_longlag <- rowSums(X_mixed)>0
    #X_trend <- setup_all_days(N)/7
    X_mu <- cbind(X_mu)
    X_M  <- cbind(X_M[,1],X_M[,2]*(X_longlag==0), X_longlag)
    p <- dim(X_mu)[2]
    p1 <- dim(X_M)[2]

    ##
    # seting up MCMC
    ##
    MH_obj <- MH_setup()
    MH_obj$sigma <- 0.05
    MH_obj$theta <- rep(0,p+p1)

    MH_obj$Lik <- loglikProbBB3
    ##
    # seting up MCMC GP
    ##
    phi.samples <- 10
    phi <- 100
    alpha.phi <- 10


    ###
    # A -matrix
    ###
    pool <- 1
    A_m  <- ceiling(N/pool)
    A_j <- rep(1:A_m,each=pool)[(A_m*pool - N + 1):(A_m*pool)]
    A   <- sparseMatrix(i=1:N,j=A_j,x=rep(1,N),  dims=c(N, A_m))
    #A[,A_m+1] <- 1


    h <- as.matrix(dist(1:N, method = "euclidean"))
    Cov <- matern.covariance(h,
                             theta[1],
                             theta[2],
                             theta[3])
    mu_Var <- t(rep(1,N))%*%solve(Cov, rep(1,N))
    Start_theta <- c(log(deaths_est+1))
    W <- solve(Cov,rep(1,N))
    mu_gp <- mean(Start_theta)
    MH_obj_GP <- MH_setup()
    MH_obj_GP$sigma <- 0.1

    MH_obj_GP$theta <- Start_theta
    MH_obj_GP$Lik <- function(x, N, phi, Sigma,A, mu) {
        prior_list <- prior_GP2(x,Sigma, mu )
        lik_list   <- lik_negbin_theta_A(N, x, phi, A)
        return(list(loglik = lik_list$loglik + prior_list$loglik,
                    grad = as.vector(lik_list$grad  + prior_list$grad)))
    }


    phis <- rep(0, MCMC_sim)

    P <- matrix(NA, ncol=N, nrow=N)
    Thetas <- matrix(NA, nrow=MCMC_sim, ncol = length(MH_obj$theta))
    theta_GP <-  matrix(NA, nrow=MCMC_sim, ncol = length(MH_obj_GP$theta))
    burnin <- ceiling(burnin_p*MCMC_sim)
    Death_est <- matrix(NA, nrow=MCMC_sim, ncol=N)
    alpha.MCMC <- rep(4, N)
    pred_set <-array(NA, dim = c(j,N_T-j,MCMC_sim))
    Alpha <- matrix(NA, ncol=N, nrow=N)
    Beta <- matrix(NA, ncol=N, nrow=N)
    mu_gps <- matrix(NA, nrow=MCMC_sim, ncol =1)
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
        MH_obj <- MALAiter(MH_obj, TRUE,
                           death.remain = data_$death.remain,
                           report.new   = data_$report.new,
                           X_mu         = X_mu,
                           X_M          = X_M)



        beta_mu    <- MH_obj$theta[1:p]
        beta_M     <- MH_obj$theta[p + (1:p1)]
        mu <- 1/(1+exp(-X_mu%*%beta_mu))
        M  <- exp(X_M%*%beta_M)
        Alpha[upper.tri(data_$report.new,diag=T)] <- M * mu
        Beta[upper.tri(data_$report.new,diag=T)]  <- M * (1-mu)

        res <- sample.deathsBB_negbin(deaths_sim,
                               deaths_est,
                               Alpha,
                               Beta,
                               Reported,
                               alpha.MCMC,
                               true.day = true.day,
                               theta = as.vector(exp(A%*%MH_obj_GP$theta)),
                               phi = phi)
        deaths_est <- res$deaths
        data_ <-newDeaths(deaths_est,Reported,maxusage.day)

        ##
        # build tau sampler
        ##
        MH_obj_GP <- MALAiter(MH_obj_GP, TRUE,
                              N = deaths_est,
                              phi = phi,
                              Sigma = Cov,
                              A = A,
                              mu = mu_gp)
        mu_gp <- rnorm(1,
                       mean = t(W)%*%MH_obj_GP$theta,
                       sd = sqrt(mu_Var))
        phi_res <- sample.phi(phi,
                              deaths_est,
                              as.vector(A%*%MH_obj_GP$theta),
                              alpha = alpha.phi)
        phi <- phi_res$phi

        if(i < burnin){
            alpha.phi <- max(1,alpha.phi + (2*(phi_res$acc/phi.samples > 0.3)  - 1))
            alpha.MCMC[res$acc/deaths_sim > 0.3] <- alpha.MCMC[res$acc/deaths_sim > 0.3] +1
            alpha.MCMC[res$acc/deaths_sim < 0.3] <- alpha.MCMC[res$acc/deaths_sim < 0.3] -1
            alpha.MCMC[alpha.MCMC<1] <- 1
        }else{
            mu_gps[i-burnin + 1,] <- mu_gp
            Thetas[i-burnin + 1,] <- MH_obj$theta
            Death_est[i-burnin + 1,] <-res$deaths
            theta_GP[i-burnin + 1,] <-  as.vector(MH_obj_GP$theta)
            phis[i-burnin + 1]  <- phi
        }
    }
    res_save <- list(Thetas = Thetas,
                     mu_gps = mu_gps,
                     Death_est = Death_est,
                     theta_GP = theta_GP,
                     true.day = true.day,
                     unique.days = unique.days,
                     phis = phis,
                     A=  A,
                     maxusage.day = maxusage.day,
                     date  = dates_report[j],
                     j = j)
    return(res_save)
}
