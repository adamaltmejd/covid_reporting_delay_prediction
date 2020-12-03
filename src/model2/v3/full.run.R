##
#


source(file.path("src","model2","v3","prob_dens.R"))
source(file.path("src","model2","v3","regression.R"))
source(file.path("src","model2","v3","model.R"))
library(foreach)
library(doParallel)
days_run <- 30
sim <- 40000
n.clust <- 3
store_data_folder <- file.path("data","tmp","model2","v3")
model_parameters <- list(sim           = sim,
                         burnin        = ceiling(0.5*sim),
                         N.days.fixed  =  3,
                         quantile      = c(0.025,0.975))

data <- readRDS(file.path("data", "processed", "processed_data.rds"))

prior_list0<- list(mu_beta        = c(0,0,0),
                   Sigma_beta     = 1/2*diag(3),
                   a_sigma        = c(3,3),
                   b_sigma        = c(5/2,5/2),
                   mu_lambda      = c(1,5),
                   Sigma_lambda   = diag(c(5,5)))

N_est_true <- apply(data$detected,1,max, na.rm=T)
N.obs <- length(data$dates)

cl <- parallel::makeCluster(n.clust)
doParallel::registerDoParallel(cl)
foreach(j = (days_run+1):(N.obs-30)) %dopar%{
    library(invgamma)
    library(mvtnorm)
    start_ = j - days_run # run the last 31 days

    report_cleaned <- report_clean(data$detected[start_:j,start_:j],data$dates[start_:j])
    new_cases <- newCases(report_cleaned)
    rownames(new_cases) <- as.character(data$dates[start_:j])
    colnames(new_cases) <- as.character(data$dates[start_:j])


    if(start_ <= days_run+ 1){
        prior_list <- prior_list0
        prior_list$n_obs <- 0
    }else{
        prior_list <- readRDS(paste(store_data_folder,"/prior_",start_-1,'.rds',sep=""))
    }

    result <- model(new_cases, model_parameters, prior_list)
    result$posteriror_list$n_obs <- prior_list$n_obs + days_run

    Npost <- cbind(data$dates[j],result$Npost, N_est_true[as.character(data$dates)%in%as.character(result$Npost$dates)])
    colnames(Npost)[c(1,7)] <- c("State",'Truth')
    saveRDS(Npost, file = paste(store_data_folder,"/Npost_",j,'.rds',sep=""))
    saveRDS(result$posteriror_list, file = paste(store_data_folder,"/prior_",j,'.rds',sep=""))
}
parallel::stopCluster(cl)


