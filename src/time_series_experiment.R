##
# simple stan code for data
#
##
library(data.table)
library(rstan)
source(file.path("src", "GPutil.R"))
source(file.path("src", "MH.R"))
sim <- 2000
rstan_options(auto_write = TRUE)
result <- readRDS(file.path("data", "processed", "processed_data.rds"))

Reported = result$detected
N <- dim(Reported)[1]

deaths_est <- apply(Reported, 1, max, na.rm=T)

death_est <- deaths_est[1:(N-15)]

N <- length(death_est)
data.list <- list(N=N,
                  Y = death_est)

plot(death_est,xlab='days',ylab='deaths')
tau <- 1000
for(i in 15 ){
    data.list_temp <- list(N=i, Y = death_est[1:i])
    stan_out = stan(file = "src/stan/RW1.stan",
                    data= data.list_temp,
                    iter = sim,
                    warmup = ceiling(0.5*sim),
                    chains =1)
    Theta_stan <- extract(stan_out)$theta
    #tau        <- extract(stan_out)$tau


    CI <- apply(exp(Theta_stan),2,quantile,c(0.1,0.5,0.9))
    lines(CI[2,],col='blue')
    lines(CI[1,],col='red')
    lines(CI[3,],col='red')
    stan_out = stan(file = "src/stan/RW1_negbin.stan",
                    data= data.list_temp,
                    iter = sim,
                    warmup = ceiling(0.5*sim),
                    chains =1)
    Theta_stan <- extract(stan_out)$theta
    #tau        <- extract(stan_out)$tau

    CI <- apply(exp(Theta_stan),2,quantile,c(0.1,0.5,0.9))
    lines(CI[2,],col='blue',lty=2)
    lines(CI[1,],col='red',lty=2)
    lines(CI[3,],col='red',lty=2)

    Start_theta <- log(data.list_temp$Y)
    ###
    # poisson part
    #
    #
    ###
    MH_obj_GP <- MH_setup()
    MH_obj_GP$sigma <- 0.01

    n_theta <- length(Start_theta)
    L <- toeplitz(c(-1,1, rep(0,n_theta-2)))
    L[lower.tri(L,diag=F)] <- 0
    L <- L[-n_theta,]
    L<-as(L, "sparseMatrix")
    MH_obj_GP$theta <- Start_theta
    MH_obj_GP$Lik <- function(x, N, L) {
        prior_list <- prior_GP(x, L)
        lik_list   <- lik_poisson(N, x)
        return(list(loglik = lik_list$loglik + prior_list$loglik,
                    grad = lik_list$grad  + prior_list$grad,
                    Hessian = lik_list$Hessian + prior_list$Hessian))
    }
    Thetas <- matrix(nrow=sim,ncol=i)
    tau_sim <- rep(0, sim)
    for(k in 1:(2*sim)){
        MH_obj_GP <- MALAiter(MH_obj_GP, TRUE,
                              N = data.list_temp$Y,
                              L = sqrt(tau) * L)
        L_theta          <- L%*%MH_obj_GP$theta
        tau       <- rgamma(1, shape=   length(L_theta)/2 + 0.01, rate = sum(L_theta^2)/2 + 0.01)
        if(k >sim){
        Thetas[k-sim,] <- as.vector(MH_obj_GP$theta)
        tau_sim[k-sim] <- tau
        }
    }
    lines(apply(exp(Thetas),2,quantile,0.5),col='green')
    lines(apply(exp(Thetas),2,quantile,0.1),col='black')
    lines(apply(exp(Thetas),2,quantile,0.9),col='black')


    ###
    # negbin part
    #
    #
    ###
    tau <- 1000
    phi <- 100
    alpha.phi <- 20
    MH_obj_GP <- MH_setup()
    MH_obj_GP$sigma <- 0.01
    n_theta <- length(Start_theta)
    L <- toeplitz(c(-1,1, rep(0,n_theta-2)))
    L[lower.tri(L,diag=F)] <- 0
    L <- L[-n_theta,]
    L<-as(L, "sparseMatrix")
    MH_obj_GP$theta <- Start_theta
    MH_obj_GP$Lik <- function(x, N, L, phi) {
        prior_list <- prior_GP(x, L)
        lik_list   <- lik_negbin_theta(N, x, phi)
        return(list(loglik = lik_list$loglik + prior_list$loglik,
                    grad = as.vector(lik_list$grad  + prior_list$grad)))
    }
    Thetas <- matrix(nrow=sim,ncol=i)
    tau_sim <- rep(0, sim)
    phi_sim <- rep(0, sim)
    for(k in 1:(2*sim)){


        MH_obj_GP <- MALAiter(MH_obj_GP, TRUE,
                              N = data.list_temp$Y,
                              L = sqrt(tau) * L,
                              phi)
        phi_res <- sample.phi(phi, data.list_temp$Y, MH_obj_GP$theta, alpha = alpha.phi)
        phi <- phi_res$phi
        L_theta          <- L%*%MH_obj_GP$theta
        tau       <- rgamma(1, shape=   length(L_theta)/2 + 0.01, rate = sum(L_theta^2)/2 + 0.01)
        if(k >sim){
            Thetas[k-sim,] <- as.vector(MH_obj_GP$theta)
            tau_sim[k-sim] <- tau
            phi_sim[k-sim] <- phi
        }else{
            alpha.phi <- max(1,alpha.phi + (2*(phi_res$acc/10 > 0.3)  - 1))
        }

    }
    lines(apply(exp(Thetas),2,quantile,0.5),col='green')
    lines(apply(exp(Thetas),2,quantile,0.1),col='black')
    lines(apply(exp(Thetas),2,quantile,0.9),col='black')
}
tau_stan <- rep(0,dim(Theta_stan)[1])
Theta_stan <- extract(stan_out)$theta


##
# testing negbin
##
test_negbin(data.list_temp$Y, extract(stan_out)$theta[1000,], extract(stan_out)$phi[1000])
X<- 1:1000
lik <- rep(0,length(X));
for(i in 1:length(X)){
    lik[i]<- MH_obj_phi$Lik(x=X[i],
                            N=data.list_temp$Y,
                            theta =  extract(stan_out)$theta[1000,])$loglik
}
plot(X,lik)

