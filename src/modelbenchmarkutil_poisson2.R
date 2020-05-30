
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
    X <- setup_data(N, maxusage.day,  dates_report, unique.days)

    ##
    # seting up MCMC
    ##
    MH_obj <- MH_setup()
    MH_obj$sigma <- 0.01
    MH_obj$theta <- rep(0,2*dim(X)[2])

    MH_obj$Lik <- loglikProbBB_param2
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
    A   <- sparseMatrix(i=1:N,j=A_j,x=rep(1,N),  dims=c(N, A_m+1))
    A[,A_m+1] <- 1


    h <- as.matrix(dist(1:N, method = "euclidean"))
    Cov_obs <- matern.covariance(h,
                             theta[1],
                             theta[2],
                             theta[3])
    Cov <- Cov
    Cov <- diag(N+1)
    Cov[1:(N),1:(N)] <- Cov_obs
    Cov[N+1,N+1]             <- theta[4]
    Start_theta <- c(log(deaths_est+1) - mean(log(deaths_est+1)),mean(log(deaths_est+1)))
    MH_obj_GP <- MH_setup()
    MH_obj_GP$sigma <- 0.1

    MH_obj_GP$theta <- Start_theta
    MH_obj_GP$Lik <- function(x, N, phi, Sigma,A) {
        prior_list <- prior_GP2(x,Sigma)
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
    p <- dim(X)[2]
    Alpha <- matrix(NA, ncol=N, nrow=N)
    Beta <- matrix(NA, ncol=N, nrow=N)
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
                           X            = X,
                           calcH= F)
        beta_1 <- MH_obj$theta[1:p]
        beta_2 <- MH_obj$theta[(p+1):(2*p)]
        mu <- 1/(1+exp(-X%*%beta_1))
        M  <- exp(X%*%beta_2)
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
                              A = A)
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
            Thetas[i-burnin + 1,] <- MH_obj$theta
            Death_est[i-burnin + 1,] <-res$deaths
            theta_GP[i-burnin + 1,] <-  as.vector(MH_obj_GP$theta)
            phis[i-burnin + 1]  <- phi
        }
    }
    res_save <- list(Thetas = Thetas,
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
