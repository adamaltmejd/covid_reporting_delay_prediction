##
#  building the benchmark prediction parameters for the full model
#  need to run Beta_GP_benchmark.R
##
library(foreach)
library(doParallel)
source(file.path("src", "util.r"))
source(file.path("src", "MH.R"))
source(file.path("src", "functions.R"))
source(file.path("src", "GPutil.R"))
nclust <- 6
MCMC_sim <- 30000
burnin_p = 0.5
deaths_sim <- 10
maxusage.day = 20 #must be less then N
unique.days  = 7
true.day = 5
start.predict.day = 26#16 # more then unique days

result <- readRDS(file.path("data", "processed", "processed_data.rds"))
Reported_T = result$detected
N_T <- dim(Reported_T)[1]

deaths_est_T <- apply(Reported_T, 1, max, na.rm=T)

data_T <- newDeaths(deaths_est_T,
                   Reported_T,
                   maxusage.day =maxusage.day)
X_T <- setup_data(N_T, maxusage.day, result$dates_report, unique.days)
predicition.list <-list()

cl <- parallel::makeCluster(nclust,setup_strategy = "sequential")
doParallel::registerDoParallel(cl)
foreach(j = start.predict.day:N_T)  %dopar% {
  library(Matrix)
  result <- readRDS(file.path("data", "processed", "processed_data.rds"))
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
  # seting up MCMC GP
  ##
  Start_theta <- log(deaths_est+1)
  tau <- 0.1
  MH_obj_GP <- MH_setup()
  MH_obj_GP$sigma <- 0.1

  MH_obj_GP$theta <- Start_theta

  MH_obj_GP$Lik <- function(x, N, tau) {
    prior_list <- prior_s1d(x, tau)
    lik_list   <- lik_poisson(N, x)
    return(list(loglik = lik_list$loglik + prior_list$loglik,
                grad = lik_list$grad  + prior_list$grad,
                Hessian = lik_list$Hessian + prior_list$Hessian,
                tau = prior_list$tau))
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
  p <- dim(X)[2]
  Alpha <- matrix(NA, ncol=N, nrow=N)
  Beta <- matrix(NA, ncol=N, nrow=N)
  for(i in 1:(MCMC_sim+burnin-1)){
    MH_obj <- MALAiter(MH_obj, TRUE,
                        death.remain = data_$death.remain,
                        report.new   = data_$report.new,
                        X            = X,
                       calcH= F)
    beta_1 <- MH_obj$theta[1:p]
    beta_2 <- MH_obj$theta[(p+1):(2*p)]
    Alpha[upper.tri(data_$report.new,diag=T)] <- exp(X%*%beta_1)
    Beta[upper.tri(data_$report.new,diag=T)]  <-  exp(X%*%beta_2)

    prior <- function(N, i){ N *MH_obj_GP$theta[i] - lgamma(N+1)}
    res <- sample.deathsBB(deaths_sim,
                           deaths_est,
                           Alpha,
                           Beta,
                           Reported,
                           alpha.MCMC,
                           true.day = true.day,
                           prior=prior)
    deaths_est <- res$deaths
    data_ <-newDeaths(deaths_est,Reported,maxusage.day)
    MH_obj_GP <- MALAiter(MH_obj_GP, TRUE,
                          N = deaths_est,
                          tau = tau)
    tau        <- MH_obj_GP$res$tau
    if(i < burnin){
      alpha.MCMC[res$acc/deaths_sim > 0.3] <- alpha.MCMC[res$acc/deaths_sim > 0.3] +1
      alpha.MCMC[res$acc/deaths_sim < 0.3] <- alpha.MCMC[res$acc/deaths_sim < 0.3] -1
      alpha.MCMC[alpha.MCMC<1] <- 1
    }else{
      Thetas[i-burnin + 1,] <- MH_obj$theta
      Death_est[i-burnin + 1,] <-res$deaths
      theta_GP[i-burnin + 1,] <-  as.vector(MH_obj_GP$theta)
    }
  }
  res_save <- list(Thetas = Thetas,
                   Death_est = Death_est,
                   theta_GP = theta_GP,
                   true.day = true.day,
                   unique.days = unique.days,
                   maxusage.day = maxusage.day)
  save(res_save,
       file = file.path("data", "simulation_results", paste0("param_", dates_report[j], ".rds")))
}

parallel::stopCluster(cl)
