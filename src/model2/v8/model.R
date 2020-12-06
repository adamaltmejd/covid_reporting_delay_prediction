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

  mu_ks    <- prior_list$mu_theta_GP
  Sigma_ks <- prior_list$Sigma_theta_GP



  mu_mu_GP    = 0
  sigma_mu_GP= 10

  mu_phi   <- prior_list$mu_phi
  sigma_phi<- prior_list$sigma_phi
  sim <- model_paramters$sim
  # start values
  sigma_time <- c(1,1)
  beta_time  <- mu_beta


  #known or fixed covariates
  #beta_fixed <- c(1)

  ###
  #
  # distribution of each death time point
  #
  ##
  Prob <- function(t){Prob_gamma(t, c(1,5))}



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
    timePoints_MH[[i]]$X <- cbind(X_time$X_first[i, index_]
                                  ,(X_time$X_t<=4)[i, index_] - X_time$X_first[i, index_],
                                  (X_time$X_t>4)[i, index_],
                                  X_time$X_na[i,index_])
    timePoints_MH[[i]]$X_fixed <- t(X_time$X_first[i, index_,drop=F])
    X_full <- rbind(X_full  , timePoints_MH[[i]]$X)
    #X_fixed <- rbind(X_fixed, timePoints_MH[[i]]$X_fixed)


    timePoints_MH[[i]]$X_sigma <- 1*cbind((X_time$X_t<=4)[i, index_],
                                          (X_time$X_t>4)[i, index_])
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

  pool <- 1
  A_m  <- ceiling(dim(new_cases)[1]/pool)
  A_j <- 1+rep(1:A_m,each=pool)[(A_m*pool - dim(new_cases)[1] + 1):(A_m*pool)]
  A   <- sparseMatrix(i=1:dim(new_cases)[1],j=A_j,  dims=c(dim(new_cases)[1], A_m+1))
  A[,1] <- 1
  MH_obj_GP <- MH_setup()
  MH_obj_GP$sigma <- 0.05

  n_theta <- dim(A)[2] + 2



  MH_obj_GP$Lik <- function(x,
                            N,
                            phi,
                            A,
                            distance_matrix,
                            mu_mu,
                            mu_sigma,
                            mu_ks,
                            Sigma_ks
                            ) {
    n <- length(x)
    kappa <- exp(x[1])
    sigma <- exp(x[2])
    Sigma_GP <- matern.covariance(distance_matrix, kappa, 1.5, sigma)
    R = chol(Sigma_GP)
    prior_ks <- dmvnorm(x[1:2], mean = mu_ks, sigma= Sigma_ks,log=T)-sum(log(diag(R)))

    prior_list <- prior_GP2(x[4:n], Sigma_GP, 0)
    prior_mu_lik <- - 0.5* (x[3]-mu_mu)^2/mu_sigma^2
    prior_mu_grad <- - (x[3]-mu_mu)/mu_sigma^2
    lik_list   <- lik_negbin_theta_A(N, x[3:n], phi, A)
    res <- list(loglik = lik_list$loglik + prior_list$loglik + prior_mu_lik + prior_ks,
                grad = c(0,0,as.vector(lik_list$grad  )))
    res$grad[3]    <- res$grad[3] +  prior_mu_grad
    res$grad[4:n]  <- res$grad[4:n]+ prior_list$grad
    return(res)
  }

  # start values
  Start_theta <- log(N+1)
  MH_obj_GP$theta <- c(0,0,mean(Start_theta), Start_theta- mean(Start_theta))
  ###
  # latent GP paramaters
  #kappa, sigma
  ###
  MH_obj_kappa_sigma <- MH_setup(Dacc_prob=0.3)
  MH_obj_kappa_sigma$sigma <- 0.05

  MH_obj_kappa_sigma$theta <- mu_ks


  d_matrix <- as.matrix(dist(1:(n_theta-3)))
  mu_ks    <- prior_list$mu_theta_GP
  Sigma_ks <- prior_list$Sigma_theta_GP

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
  theta_vec        <- matrix(0, nrow=sim,   ncol=n_theta-3)
  phi_vec          <- matrix(0, nrow = sim, ncol = 1)
  Beta_vec         <- matrix(0, nrow=sim,   ncol=length(mu_beta))
  sigma_vec        <- matrix(0, nrow=sim,   ncol=2)
  ProbMatrix_vec   <- matrix(0, nrow=sim,   ncol= N_0*N_0)
  sigma_theta_vec  <- matrix(0, nrow=sim,  ncol = 1)
  kappa_vec        <- matrix(0, nrow=sim,  ncol = 1)
  mu_GP_vec        <- matrix(0, nrow=sim,  ncol = 1)
  for(iter in 1:sim){

    # simulating time pointsm N and hyperparameters
    report_effort <- c()
    ProbMatrix_new <- matrix(NA, N_0, N_0)
    lambda <- exp(A%*%MH_obj_GP$theta[3:n_theta]) #prior on N from GP
    for(i in 1:dim(new_cases)[1]){
      timePoints_MH[[i]] <- MHiter(timePoints_MH[[i]],calcLik = T,
                                   timePoints_MH[[i]]$n.obs,
                                   N = N[i],
                                   timePoints_MH[[i]]$X%*%beta_time,# + timePoints_MH[[i]]$X_fixed%*%beta_fixed,
                                   timePoints_MH[[i]]$X_sigma%*%sigma_time,
                                   Prob)
      report_effort <- c(report_effort, timePoints_MH[[i]]$theta)
      ProbMatrix_new[i, timePoints_MH[[i]]$index_non_na] <- timePoints_MH[[i]]$theta
      if(i > model_paramters$N.days.fixed){
        # sampling N

        Nstar <- sample((N[i]-alpha.MCMC[i]):(N[i]+alpha.MCMC[i] ), 1)
        if(Nstar >= sum(timePoints_MH[[i]]$n.obs)){
          lik      <- density_t((timePoints_MH[[i]]$theta), timePoints_MH[[i]]$n.obs, N[i], Prob)   + dnegbin(N[i], lambda[i],phi)
          lik_star <-  density_t((timePoints_MH[[i]]$theta), timePoints_MH[[i]]$n.obs, Nstar, Prob) + dnegbin(Nstar, lambda[i],phi)
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
                          distance_matrix = d_matrix,
                          mu_mu    = mu_mu_GP,
                          mu_sigma = sigma_mu_GP,
                          mu_ks =mu_ks,
                          Sigma_ks = Sigma_ks)
    kappa <- exp(MH_obj_GP$theta)[1]
    sigma2_theta <- exp(MH_obj_GP$theta)[2]

    phi_res <- sample.phi.prior(phi,
                                N,
                                c(mu_phi, sigma_phi),
                                as.vector(A%*%MH_obj_GP$theta[3:n_theta]),
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
      sigma_time    <- sqrt(rinvgamma(2, shape= a_hat, rate=b_hat))


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
    theta_vec[iter, ]   <- as.vector(A%*%as.vector(MH_obj_GP$theta[3:n_theta]))
    phi_vec[iter, ]     <- phi
    sigma_theta_vec[iter,]   <- sqrt(sigma2_theta)
    kappa_vec[iter, ]         <- kappa
    mu_GP_vec[iter, ]    <-MH_obj_GP$theta[3]

  }
  Beta_vec        <- Beta_vec[model_paramters$burnin:model_paramters$sim, ]
  sigma_vec       <- sigma_vec[model_paramters$burnin:model_paramters$sim,]
  N_vec           <- N_vec[model_paramters$burnin:model_paramters$sim,]
  ProbMatrix_vec  <- ProbMatrix_vec[model_paramters$burnin:model_paramters$sim,]
  phi_vec         <- phi_vec[model_paramters$burnin:model_paramters$sim,]
  theta_vec       <- theta_vec[model_paramters$burnin:model_paramters$sim,]
  sigma_theta_vec <- sigma_theta_vec[model_paramters$burnin:model_paramters$sim,]
  kappa_vec       <- kappa_vec[model_paramters$burnin:model_paramters$sim,]
  mu_GP_vec       <- mu_GP_vec[model_paramters$burnin:model_paramters$sim,]
  post_sigma_1 <- ml_inversegamma(sigma_vec[,1]^2)
  post_sigma_2 <- ml_inversegamma(sigma_vec[,2]^2)


  posterior_list <- list(mu_beta      = colMeans(Beta_vec),
                         Sigma_beta   = cov(Beta_vec),
                         a_sigma      = c(post_sigma_1[1],post_sigma_2[1]),
                         b_sigma      = c(post_sigma_1[2],post_sigma_2[2]),
                         mu_phi       = mean(log(phi_vec)),
                         sigma_phi       = sd(log(phi_vec)),
                         mu_theta_GP      = colMeans(log(cbind(kappa_vec,sigma_theta_vec^2))),
                         Sigma_theta_GP   = cov(log(cbind(kappa_vec,sigma_theta_vec^2))))


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
                                   sigma_theta = sigma_theta_vec,
                                   kappa       = kappa_vec,
                                   mu_GP       = mu_GP_vec)
  result$X_full <- X_full
  return(result)
}





