library(tidyr)
source(file.path("src", "util","util.r"))
source(file.path("src", "util","MH.R"))
source(file.path("src", "util","functions.R"))
source(file.path("src", "util","GPutil.R"))
source(file.path("src","util","MLbeta.R"))
library(numDeriv)
library(invgamma)
library(mvtnorm)


#'
#' model for estimating number of deaths
#' Each death is given an effort and
#' each days one works a random effort on each death once work>effort
#' the death is reported
#' No prior on number of deaths is given
#'
#' @param data -  (TxT ) of new_cases
#' @param model_parameters -
#'                          [sim]     number of MCMC samples
#'                          [burnin]  discaredd number of MCMC samples
#'                          [N.days.fixed] - the first N.days.fixed number of death is the reported cases
#'                          [quantile]     - lower uper confidence intervall for deaths
#' @param prior_list       - list of priors for the parameters
#'                         - [mu_beta]    mean     fixed effect parameters
#'                         - [Sigma_beta] variance fixed effect parameters
#'                         - [a_sigma]    variance parameter for the first seven days
#'                         - [b_sigma]    variance parameter after seven days
model <- function(new_cases, model_paramters, prior_list, startvalue_list=NULL){

  N_0 <- dim(new_cases)[1]
  dates <- colnames(new_cases)
  mu_beta <- prior_list$mu_beta
  Q_beta <- solve(prior_list$Sigma_beta)
  #sigma prior
  a_sigma <- prior_list$a_sigma
  b_sigma <- prior_list$b_sigma

  theta_mu    <- prior_list$theta_mu
  theta_Sigma <- prior_list$theta_Sigma




  mu_phi   <- prior_list$mu_phi
  sigma_phi<- prior_list$sigma_phi
  sim <- model_paramters$sim
  # start values
  sigma_time <- b_sigma/(a_sigma-1)
  beta_time  <- mu_beta


  #known or fixed covariates
  #beta_fixed <- c(1)

  ###
  #
  # distribution of each death time point
  #
  ##
  Prob <- function(t){Prob_gamma(t, c(1,10))}



  ##
  # build models for time point distribution,
  # and sampling N
  ###
  X_time <- build_X_time(new_cases)
  new_cases[X_time$X_diag==1] <- NA # remove zero first observation

  dates_keep <- rowSums(is.na(new_cases))<dim(new_cases)[2]
  new_cases <- new_cases[dates_keep,]



  if(is.null(startvalue_list))
    N <-apply(new_cases,1,sum, na.rm=T)
  X_full <- c()
  X_fixed <- c()
  X_sigma <- c()
  timePoints_MH <- list()
  n.obs.full <- c()
  max_t <- 0 #largest number of ts
  for(i in 1:dim(new_cases)[1]){
    cases_i <- new_cases[i,]
    index_  <- is.na(cases_i)==F
    max_t <- max(c(max_t, sum(index_)))

    timePoints_MH[[i]] <- MH_setup(Dacc_prob=0.3)
    timePoints_MH[[i]]$sigma <- 0.05
    timePoints_MH[[i]]$theta <- rep(1, sum(index_))
    timePoints_MH[[i]]$n.obs <- cases_i[index_]
    timePoints_MH[[i]]$index_non_na <-      index_
    timePoints_MH[[i]]$Lik <- function(x,
                                       n.obs,
                                       N,
                                       mu_t,
                                       sigma_t,
                                       Prob) {
      if(min(x)<0)
        return(list(loglik =-Inf))
      log_lik   <- density_t(x, n.obs, N, Prob)
      lik_prior  <- prior_t(x,mu_t, sigma_t)
      return(list(loglik = log_lik + lik_prior))
    }
    n.obs.full <- c(n.obs.full,timePoints_MH[[i]]$n.obs)


    #setting upp covariates for priro
    timePoints_MH[[i]]$X <- cbind(X_time$first[i, index_],
                                  X_time$second[i, index_] +
                                      X_time$fourth[i, index_]+
                                  X_time$third[i, index_],
                                  X_time$rest[i, index_],
                                  X_time$Tuesday[i,index_])
    X_full <- rbind(X_full  , timePoints_MH[[i]]$X)
    #X_fixed <- rbind(X_fixed, timePoints_MH[[i]]$X_fixed)


    timePoints_MH[[i]]$X_sigma <- 1*cbind(X_time$first[i, index_],
                                           X_time$second[i, index_] +
                                               X_time$fourth[i, index_]+
                                               X_time$third[i, index_],
                                           X_time$rest[i, index_])
    X_sigma <- rbind(X_sigma,
                     timePoints_MH[[i]]$X_sigma)
    if(i > dim(new_cases)[1]-14){
      N[i] <- max(ceiling(mean(cases_i[index_][1:min(2,sum(index_))])/0.3), N[i])
      }
  }
  # for sampling adaptive N sampling
  alpha.MCMC <- rep(10, dim(new_cases)[1])
  acc_N      <- rep(0, dim(new_cases)[1])


  XXt <- t(X_full)%*%X_full
  Q_post     <- XXt +  Q_beta
  Sigma_post <- solve(Q_post)
  L_post <- t(chol(Sigma_post))


  ####
  #
  # GP prior
  #
  ####

  A   <- matrix(0, nrow= dim(new_cases)[1], ncol = 2)
  A[,1] <- 1
  A[,2] <- 1:dim(new_cases)[1]
  MH_obj_GP <- MH_setup()
  MH_obj_GP$sigma <- 0.005
    n_theta <- 2

  MH_obj_GP$Lik <- function(x,
                            N,
                            phi,
                            A,
                            mu_theta,
                            sigma_theta
                            ) {


    Sigma_x    <-   solve(sigma_theta,x-mu_theta)
    lik_prior  <-  -0.5  * t(x-mu_theta)%*%Sigma_x
    grad_prior <-  -Sigma_x
    lik_list   <- lik_negbin_theta_A(N, x, phi, A)
    res <- list(loglik = lik_list$loglik + lik_prior,
                grad = lik_list$grad + grad_prior )
    return(res)
  }

  Start_theta <- log(N+1)
  MH_obj_GP$theta <- c(mean(Start_theta), 0)


  # Neg bin part
  phi <- 10 # start value
  alpha.phi <- 10 #AMCMC param
  phi.samples <- 10 #AMCMC param
  phi.samples.total <- 0
  acc_phi <- 0
  ####
  # GP prior done
  ####

  N_vec <- matrix(0,nrow=sim, ncol=length(N))
  colnames(N_vec) <- rownames(new_cases)
  theta_vec        <- matrix(0, nrow=sim,   ncol=n_theta)
  phi_vec          <- matrix(0, nrow = sim, ncol = 1)
  Beta_vec         <- matrix(0, nrow=sim,   ncol=length(mu_beta))
  sigma_vec        <- matrix(0, nrow=sim,   ncol=length(a_sigma))
  ProbMatrix_vec   <- matrix(0, nrow=sim,   ncol= N_0*N_0)
  t_vec            <- matrix(NA, nrow=sim, ncol= dim(new_cases)[1] * max_t)
  for(iter in 1:sim){

    # simulating time pointsm N and hyperparameters
    report_effort <- c()
    ProbMatrix_new <- matrix(NA, N_0, N_0)
    lambda <- exp(A%*%MH_obj_GP$theta) #prior on N from GP
    for(i in 1:dim(new_cases)[1]){
      timePoints_MH[[i]] <- MHiter(timePoints_MH[[i]],calcLik = T,
                                   timePoints_MH[[i]]$n.obs,
                                   N = N[i],
                                   timePoints_MH[[i]]$X%*%beta_time,# + timePoints_MH[[i]]$X_fixed%*%beta_fixed,
                                   timePoints_MH[[i]]$X_sigma%*%sigma_time,
                                   Prob)
      report_effort <- c(report_effort, timePoints_MH[[i]]$theta)
      t_vec[iter, max_t*(i-1)+ (1:sum(timePoints_MH[[i]]$index_non_na))] <- (timePoints_MH[[i]]$theta - timePoints_MH[[i]]$X%*%beta_time)
      ProbMatrix_new[i, timePoints_MH[[i]]$index_non_na] <- Prob(timePoints_MH[[i]]$theta)[1:(length(timePoints_MH[[i]]$theta))]
      if(i > model_paramters$N.days.fixed){
        # sampling N

        Nstar <- sample((N[i]-alpha.MCMC[i]):(N[i]+alpha.MCMC[i] ), 1)
        if(Nstar >= sum(timePoints_MH[[i]]$n.obs)){
          lik      <- density_t((timePoints_MH[[i]]$theta), timePoints_MH[[i]]$n.obs, N[i], Prob)   #+ dnegbin(N[i], lambda[i],phi)
          lik_star <-  density_t((timePoints_MH[[i]]$theta), timePoints_MH[[i]]$n.obs, Nstar, Prob) #+ dnegbin(Nstar, lambda[i],phi)
          if(log(runif(1)) < lik_star-lik ){
            N[i] <- Nstar
            acc_N[i] <- acc_N[i] + 1
          }
        }
      }
    }


    MH_obj_GP <- MALAiter(MH_obj_GP, TRUE,
                          N = N,
                          phi = phi,
                          A = A,
                          mu_theta = theta_mu,
                          sigma_theta = theta_Sigma)


    phi_res <- sample.phi.prior(phi,
                                N,
                                c(mu_phi, sigma_phi),
                                as.vector(A%*%MH_obj_GP$theta),
                                alpha = alpha.phi,
                                samples = phi.samples)
    phi <- phi_res$phi
    acc_phi <- acc_phi + phi_res$acc
    phi.samples.total <- phi.samples.total + phi.samples


    ##
    # sampling latent parameters
    #
    ###
    if(iter >100){
      Q_obs      <- diag(1./as.vector(X_sigma%*%sigma_time)^2)
      Q_post     <- t(X_full)%*%Q_obs%*%(X_full) +  Q_beta
      Sigma_post <- solve(Q_post)

      L_post <- t(chol(Sigma_post))
      #report_effort <-  report_effort -  X_fixed%*%beta_fixed
      mu_hat <- Sigma_post%*%(t(X_full)%*%(Q_obs%*%report_effort) + Q_beta%*%mu_beta)

      beta_time     <- mu_hat +   ( L_post%*%rnorm(length(beta_time)))

      a_hat  <- a_sigma + colSums(X_sigma)/2
      report_effort = report_effort - X_full%*%beta_time
      b_hat  <- b_sigma + 0.5  *t(X_sigma)%*%(report_effort)^2
      sigma_time    <- sqrt(rinvgamma(length(a_hat), shape= a_hat, rate=b_hat))


    }



    if(iter%%50==0 &  iter < model_parameters$burnin){
      alpha.phi <- max(1,alpha.phi + (2*(acc_phi/phi.samples.total > 0.3)  - 1))
      alpha.MCMC[acc_N/50 > 0.3] <- alpha.MCMC[acc_N/50 > 0.3] +1
      alpha.MCMC[acc_N/50 < 0.3] <- alpha.MCMC[acc_N/50 < 0.3] -1
      alpha.MCMC[alpha.MCMC<1] <- 1
      acc_N             <-acc_N* 0
      phi.samples.total <- 0
      acc_phi           <- 0
    }

    Beta_vec[iter, ] <- beta_time #iXXt%*%t(X_full)%*%log_t
    sigma_vec[iter,] <- sigma_time
    ProbMatrix_vec[iter,] <- ProbMatrix_new
    N_vec[iter,] <- N
    theta_vec[iter, ]   <- MH_obj_GP$theta
    phi_vec[iter, ]     <- phi

  }
  Beta_vec        <- Beta_vec[model_paramters$burnin:model_paramters$sim, ]
  sigma_vec       <- sigma_vec[model_paramters$burnin:model_paramters$sim,]
  N_vec           <- N_vec[model_paramters$burnin:model_paramters$sim,]
  ProbMatrix_vec  <- ProbMatrix_vec[model_paramters$burnin:model_paramters$sim,]
  phi_vec         <- phi_vec[model_paramters$burnin:model_paramters$sim,]
  theta_vec       <- theta_vec[model_paramters$burnin:model_paramters$sim,]
  t_vec           <- t_vec[model_paramters$burnin:model_paramters$sim,]
  a_vec <- c()
  b_vec <- c()
  for(i in 1:dim(sigma_vec)){
      post_sigma_1 <- ml_inversegamma(sigma_vec[,1]^2)
      a_vec <- c(a_vec,post_sigma_1[1])
      b_vec <- c(b_vec,post_sigma_1[2])
  }

  post_sigma_2 <- ml_inversegamma(sigma_vec[,2]^2)


  posterior_list <- list(mu_beta      = colMeans(Beta_vec),
                         Sigma_beta   = cov(Beta_vec),
                         a_sigma      = a_vec,
                         b_sigma      = b_vec,
                         mu_phi       = mean(log(phi_vec)),
                         sigma_phi       = sd(log(phi_vec)),
                         theta_mu      = colMeans(theta_vec),
                         theta_Sigma   = cov(theta_vec))


  result <- list()
  CI <- t(apply(N_vec,2 , quantile, prob=model_parameters$quantile))
  result$Npost <- data.frame(dates = rownames(new_cases),
                             mean  = colMeans(N_vec),
                             median  =apply(N_vec,2 , quantile, prob=0.5),
                             lCI   = CI[,1],
                             uCI   = CI[,2])

  result$posteriror_list <- posterior_list
  result$posteriror_sample <- list(Beta        = Beta_vec,
                                   sigma      = sigma_vec,
                                   N          = N_vec,
                                   theta      = theta_vec,
                                   phi        = phi_vec,
                                   ProbMatrix = ProbMatrix_vec,
                                   t          = t_vec,
                                   max_t      = max_t)
  result$X_full <- X_full
  return(result)
}





