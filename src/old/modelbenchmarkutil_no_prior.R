
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
##
benchmark_BetaGP_j <- function(j,
                               result,
                               true.day,
                               maxusage.day,
                               unique.days,
                               MCMC_sim,
                               burnin_p,
                               deaths_sim){
    require(Matrix)
    prior_no <- function(x,i) {dnbinom(x,size=4,prob=0.035,log=T)}
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
    MH_obj$sigma <- 0.1
    MH_obj$theta <- rep(0,2*dim(X)[2])


    MH_obj$Lik <- loglikProbBB


    ##
    # mcmc loop
    ##
    P <- matrix(NA, ncol=N, nrow=N)
    Thetas <- matrix(NA, nrow=MCMC_sim, ncol = length(MH_obj$theta))
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
        Alpha[upper.tri(data_$report.new,diag=T)] <- exp(X%*%beta_1)
        Beta[upper.tri(data_$report.new,diag=T)]  <-  exp(X%*%beta_2)

        res <- sample.deathsBB(deaths_sim,
                               deaths_est,
                               Alpha,
                               Beta,
                               Reported,
                               alpha.MCMC,
                               true.day = true.day,
                               use.prior=T,
                               prior=prior_no)
        deaths_est <- res$deaths
        data_ <-newDeaths(deaths_est,Reported,maxusage.day)

        ##
        # build tau sampler
        ##

        if(i < burnin){
            alpha.MCMC[res$acc/deaths_sim > 0.3] <- alpha.MCMC[res$acc/deaths_sim > 0.3] +1
            alpha.MCMC[res$acc/deaths_sim < 0.3] <- alpha.MCMC[res$acc/deaths_sim < 0.3] -1
            alpha.MCMC[alpha.MCMC<1] <- 1
        }else{
            Thetas[i-burnin + 1,] <- MH_obj$theta
            Death_est[i-burnin + 1,] <-res$deaths
        }
    }
    res_save <- list(Thetas = Thetas,
                     Death_est = Death_est,
                     true.day = true.day,
                     unique.days = unique.days,
                     maxusage.day = maxusage.day)
    return(res_save)
}
