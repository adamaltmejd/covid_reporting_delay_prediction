
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
#' prior        - (2 x 1) [1] 0: rw1, 1: rw2,
#'                        [2] 0: no noise, 1: noise
##
benchmark_BinomGP_j <- function(j,
                               result,
                               true.day,
                               maxusage.day,
                               unique.days,
                               MCMC_sim,
                               burnin_p,
                               deaths_sim,
                               prior = c(0,0)){
    require(Matrix)
    N_T <- dim(result$detected)[1]
    deaths_est_T <- apply(result$detected, 1, max, na.rm=T)
    Reported <- result$detected[1:j,1:j]
    dates_report <- result$dates_report[1:j]
    N <- dim(Reported)[1]
    deaths_est <- apply(Reported, 1, max, na.rm=T)
    deaths_est[1:true.day] = deaths_est_T[1:true.day] #days with known deaths


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
    MH_obj$theta <- rep(0,dim(X)[2])


    MH_obj$Lik <- loglikProb


    ##
    # seting up MCMC GP
    ##

    phi.samples <- 10
    phi <- 100
    alpha.phi <- 10
    Start_theta <- log(deaths_est+1)
    tau <- 0.1
    MH_obj_GP <- MH_setup()
    MH_obj_GP$sigma <- 0.1

    MH_obj_GP$theta <- Start_theta
    if(prior[1] == 0){
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
    Gp.prior <- matrix(NA, nrow=MCMC_sim, ncol=1)
    if(prior[2] == 1){
        #MH sampler for hyperpparameters
        MH_obj_tau <- MH_setup(Dacc_prob=0.3)
        MH_obj_tau$sigma <- 0.1
        MH_obj_tau$Lik <- tau_sigma_prior
        MH_obj_tau$theta <- c(tau, 1)
        Gp.prior <- matrix(NA, nrow=MCMC_sim, ncol=2)
    }
    MH_obj_GP$Lik <- function(x, N, L, phi) {
        prior_list <- prior_GP(x, L)
        lik_list   <- lik_negbin_theta(N, x, phi)
        return(list(loglik = lik_list$loglik + prior_list$loglik,
                    grad = as.vector(lik_list$grad  + prior_list$grad)))
    }

    ##
    # mcmc loop
    ##
    P <- matrix(NA, ncol=N, nrow=N)
    Thetas <- matrix(NA, nrow=MCMC_sim, ncol = length(MH_obj$theta))
    theta_GP <-  matrix(NA, nrow=MCMC_sim, ncol = length(MH_obj_GP$theta))
    burnin <- ceiling(burnin_p*MCMC_sim)
    Death_est <- matrix(NA, nrow=MCMC_sim, ncol=N)
    alpha.MCMC <- rep(4, N)
    pred_set <-array(NA, dim = c(j,N_T-j,MCMC_sim))
    phis <- rep(0, MCMC_sim)
    p <- dim(X)[2]
    P_avg <- matrix(0, ncol=N, nrow=N)
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
                           X            = X)
        P[upper.tri(data_$report.new,diag=T)] <- as.vector(1/(1+ exp(-X%*%MH_obj$theta)))
        prior_N <- function(N, i){ dnegbin(N, exp(MH_obj_GP$theta[i]),phi) }
        res <- sample.deaths(deaths_sim,
                               deaths_est,
                               P,
                               Reported,
                               alpha.MCMC,
                               true.day = true.day,
                               use.prior=T,
                               prior=prior_N)
        deaths_est <- res$deaths
        data_ <-newDeaths(deaths_est,Reported,maxusage.day)

        ##
        # build tau sampler
        ##

        if(prior[2] == 0){
            L_theta          <- L%*%MH_obj_GP$theta

            tau   <- rgamma(1,
                            shape=length(L_theta)/2 + 0.01,
                            scale = sum(L_theta^2)/2 + 0.01)
            L_in <- sqrt(tau) * L
        }else{
            L_theta          <- (E_Q%*%MH_obj_GP$theta)
            MH_obj_tau <- MHiter(MH_obj_tau,
                                 calcLik=TRUE,
                                 L_theta=  L_theta,
                                 Q_E = D_Q)
            L_in <- diag(sqrt(1/(1/(MH_obj_tau$theta[1]*D_Q) + 1/MH_obj_tau$theta[2])))%*%E_Q
        }

        MH_obj_GP <- MALAiter(MH_obj_GP, TRUE,
                              N = deaths_est,
                              L = L_in,
                              phi = phi)
        phi_res <- sample.phi(phi,
                              deaths_est,
                              MH_obj_GP$theta,
                              alpha = alpha.phi)
        phi <- phi_res$phi
        if(i < burnin){
            alpha.phi <- max(1,alpha.phi + (2*(phi_res$acc/phi.samples > 0.3)  - 1))
            alpha.MCMC[res$acc/deaths_sim > 0.3] <- alpha.MCMC[res$acc/deaths_sim > 0.3] +1
            alpha.MCMC[res$acc/deaths_sim < 0.3] <- alpha.MCMC[res$acc/deaths_sim < 0.3] -1
            alpha.MCMC[alpha.MCMC<1] <- 1
        }else{
            P_avg <- P_avg + P/MCMC_sim
            Thetas[i-burnin + 1,] <- MH_obj$theta
            Death_est[i-burnin + 1,] <-res$deaths
            theta_GP[i-burnin + 1,] <-  as.vector(MH_obj_GP$theta)
            if(prior[2]==0){
                Gp.prior[i-burnin + 1]  <- tau
            }else{
                Gp.prior[i-burnin + 1,]  <- MH_obj_tau$theta

            }
            phis[i-burnin + 1] <- phi
        }
    }
    res_save <- list(Thetas = Thetas,
                     Death_est = Death_est,
                     theta_GP = theta_GP,
                     true.day = true.day,
                     unique.days = unique.days,
                     maxusage.day = maxusage.day,
                     GP.prior      = Gp.prior,
                     phis  =  phis,
                     P_avg = P_avg )
    return(res_save)
}
