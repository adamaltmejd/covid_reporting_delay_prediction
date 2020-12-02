library(tidyr)
source(file.path("src", "util","util.r"))
source(file.path("src", "util","MH.R"))
source(file.path("src", "util","functions.R"))
source(file.path("src", "util","GPutil.R"))
source(file.path("src","util","MLbeta.R"))
library(numDeriv)
library(invgamma)
library(mvtnorm)

# For version 2
# adapt N samples


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

  sim <- model_paramters$sim
  # start values
  sigma_time <- c(1,1)
  beta_time  <- mu_beta


  #known or fixed covariates
  beta_fixed <- c(1)

  ###
  #
  # distribution of each death time point
  #
  ##
  effort_MH <- MH_setup(Dacc_prob=0.3)
  Gamma_theta <- c(5)
  effort_MH$sigma <- 0.1
  effort_MH$theta <- Gamma_theta
  effort_MH$Lik <- function(x,
                            Time_objs,
                            N,
                            mu,
                            Sigma) {
    x <- as.vector(x)
    if(x<0)
      return(list(loglik =-Inf))
    log_lik <- 0
    Prob <- function(t){Prob_gamma(t, c(1,x))}
    for(i in 1:length(Time_objs)){
      log_lik   <- log_lik +  density_t(Time_objs[[i]]$theta, Time_objs[[i]]$n.obs, N[i], Prob)
    }
    lik_prior  <- dnorm(x, mu, Sigma, log=T)
    return(list(loglik = log_lik + lik_prior))
  }
  Prob <- function(t){Prob_gamma(t, c(1,effort_MH$theta))}



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
  for(i in 1:dim(new_cases)[1]){
    cases_i <- new_cases[i,]
    index_  <- is.na(cases_i)==F
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
    timePoints_MH[[i]]$X <- cbind((X_time$X_t<=4)[i, index_] - X_time$X_first[i, index_],
                                  (X_time$X_t>4)[i, index_],
                                  X_time$X_na[i,index_])
    timePoints_MH[[i]]$X_fixed <- t(X_time$X_first[i, index_,drop=F])
    X_full <- rbind(X_full  , timePoints_MH[[i]]$X)
    X_fixed <- rbind(X_fixed, timePoints_MH[[i]]$X_fixed)


    timePoints_MH[[i]]$X_sigma <- 1*cbind((X_time$X_t<=4)[i, index_],
                                          (X_time$X_t>4)[i, index_])
    X_sigma <- rbind(X_sigma,
                     timePoints_MH[[i]]$X_sigma)

  }
  # for sampling adaptive N sampling
  alpha.MCMC <- rep(10, dim(new_cases)[1])
  acc_N      <- rep(0, dim(new_cases)[1])


  XXt <- t(X_full)%*%X_full
  Q_post     <- XXt +  Q_beta
  Sigma_post <- solve(Q_post)
  L_post <- t(chol(Sigma_post))



  N_vec <- matrix(0,nrow=sim, ncol=length(N))
  colnames(N_vec) <- rownames(new_cases)
  Beta_vec       <- matrix(0,nrow=sim,    ncol=length(mu_beta))
  sigma_vec      <- matrix(0,nrow=sim,   ncol=2)
  lambda_vec     <- matrix(0, nrow=sim, ncol=1)
  ProbMatrix_vec <- matrix(0, nrow=sim, ncol= N_0*N_0)
  for(iter in 1:sim){

    # simulating time pointsm N and hyperparameters
    report_effort <- c()
    ProbMatrix_new <- matrix(NA, N_0, N_0)
    for(i in 1:dim(new_cases)[1]){
      timePoints_MH[[i]] <- MHiter(timePoints_MH[[i]],calcLik = T,
                                   timePoints_MH[[i]]$n.obs,
                                   N = N[i],
                                   timePoints_MH[[i]]$X%*%beta_time + timePoints_MH[[i]]$X_fixed%*%beta_fixed, #add plus one for first occurence
                                   timePoints_MH[[i]]$X_sigma%*%sigma_time,
                                   Prob)
      report_effort <- c(report_effort, timePoints_MH[[i]]$theta)
      ProbMatrix_new[i, timePoints_MH[[i]]$index_non_na] <- timePoints_MH[[i]]$theta
      #ProbMatrix_new[i, timePoints_MH[[i]]$index_non_na] <- Prob(timePoints_MH[[i]]$theta)[1:length(timePoints_MH[[i]]$theta)]
      if(i > model_paramters$N.days.fixed){
        # sampling N

        Nstar <- sample((N[i]-alpha.MCMC[i]):(N[i]+alpha.MCMC[i] ), 1)
        if(Nstar >= sum(timePoints_MH[[i]]$n.obs)){
          lik      <- density_t((timePoints_MH[[i]]$theta), timePoints_MH[[i]]$n.obs, N[i], Prob) + dnegbin(N[i], 50, 2)
          lik_star <-  density_t((timePoints_MH[[i]]$theta), timePoints_MH[[i]]$n.obs, Nstar, Prob)  + dnegbin(Nstar, 50, 2)
          if(log(runif(1)) < lik_star-lik ){
            N[i] <- Nstar
            acc_N[i] <- acc_N[i] + 1
          }
        }
      }
    }
    effort_MH <- MHiter(effort_MH,
                        calcLik = T,
                        Time_objs = timePoints_MH,
                        N = N,
                        mu = prior_list$mu_lambda,
                        Sigma = prior_list$Sigma_lambda
                        )
    Prob <- function(t){Prob_gamma(t, c(1,effort_MH$theta))}
    N_vec[iter,] <- N
    ##
    # sampling latent parameters
    #
    ###
    if(iter >100){
      Q_obs      <- diag(1./as.vector(X_sigma%*%sigma_time)^2)
      Q_post     <- t(X_full)%*%Q_obs%*%(X_full) +  Q_beta
      Sigma_post <- solve(Q_post)

      L_post <- t(chol(Sigma_post))
      report_effort <-  report_effort -  X_fixed%*%beta_fixed
      mu_hat <- Sigma_post%*%(t(X_full)%*%(Q_obs%*%report_effort) + Q_beta%*%mu_beta)

      beta_time     <- mu_hat +   ( L_post%*%rnorm(length(beta_time)))

      a_hat  <- a_sigma + colSums(X_sigma)/2
      report_effort = report_effort - X_full%*%beta_time
      b_hat  <- b_sigma + 0.5  *t(X_sigma)%*%(report_effort)^2
      sigma_time    <- sqrt(rinvgamma(2, shape= a_hat, rate=b_hat))

    }


    if(iter%%50==0&  iter < model_parameters$burnin){
      alpha.MCMC[acc_N/50 > 0.3] <- alpha.MCMC[acc_N/50 > 0.3] +1
      alpha.MCMC[acc_N/50 < 0.3] <- alpha.MCMC[acc_N/50 < 0.3] -1
      alpha.MCMC[alpha.MCMC<1] <- 1
      acc_N <-acc_N* 0
    }

    Beta_vec[iter, ] <- beta_time #iXXt%*%t(X_full)%*%log_t
    sigma_vec[iter,] <- sigma_time
    lambda_vec[iter, ] <- effort_MH$theta
    ProbMatrix_vec[iter,] <- ProbMatrix_new
  }
  Beta_vec   <- Beta_vec[model_paramters$burnin:model_paramters$sim, ]
  sigma_vec  <- sigma_vec[model_paramters$burnin:model_paramters$sim,]
  N_vec      <- N_vec[model_paramters$burnin:model_paramters$sim,]
  lambda_vec <- lambda_vec[model_paramters$burnin:model_paramters$sim]
  ProbMatrix_vec <- ProbMatrix_vec[model_paramters$burnin:model_paramters$sim,]
  post_sigma_1 <- ml_inversegamma(sigma_vec[,1])
  post_sigma_2 <- ml_inversegamma(sigma_vec[,2])
  posterior_list <- list(mu_beta      = colMeans(Beta_vec),
                         Sigma_beta   = cov(Beta_vec),
                         a_sigma      = c(post_sigma_1[1],post_sigma_2[1]),
                         b_sigma      = c(post_sigma_1[2],post_sigma_2[2]),
                         mu_lambda    = mean(lambda_vec),
                         Sigma_lambda = sd(lambda_vec))

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
                                   lambda     = lambda_vec ,
                                   N          = N_vec,
                                   ProbMatrix = ProbMatrix_vec)
  return(result)
}





