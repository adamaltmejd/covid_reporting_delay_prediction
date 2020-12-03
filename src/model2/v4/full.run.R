##
#


source(file.path("src","model2","v4","prob_dens.R"))
source(file.path("src","model2","v4","regression.R"))
source(file.path("src","model2","v4","model.R"))
library(foreach)
library(doParallel)
days_run <- 30
sim <- 40000
n.clust <- 3
store_data_folder <- file.path("data","tmp","model2","v4")
model_parameters <- list(sim           = sim,
                         burnin        = ceiling(0.5*sim),
                         N.days.fixed  =  3,
                         quantile      = c(0.025,0.975))

data <- readRDS(file.path("data", "processed", "processed_data.rds"))

prior_list <- list(mu_beta        = c(0,0,0),
                   Sigma_beta     = 1/2*diag(3),
                   a_sigma        = c(3,3),
                   b_sigma        = c(5/2,5/2),
                   mu_lambda      = 5,
                   Sigma_lambda   = 5,
                   mu_GP          = 0,
                   sigma_GP       = 10,
                   mu_phi         = 1,
                   sigma_phi      = 1,
                   a_sigma_theta  = 4,
                   b_sigma_theta  = 0.1*(4-1),
                   sigma2_theta_max = 0.4806756^2)

N_est_true <- apply(data$detected,1,max, na.rm=T)
N.obs <- length(data$dates)

cl <- parallel::makeCluster(n.clust)
doParallel::registerDoParallel(cl)
foreach(j = (days_run+1):(N.obs-30)) %dopar%{
    library(invgamma)
    start_ = j - days_run # run the last 31 days

    if(start_ <= days_run+ 1){
        prior_list <- prior_list0
        prior_list$n_obs <- 0
        start_ = 1
    }else{
        prior_list <- readRDS(paste(store_data_folder,"/prior_",start_-1,'.rds',sep=""))
    }

    report_cleaned <- report_clean(data$detected[start_:j,start_:j],data$dates[start_:j])
    new_cases <- newCases(report_cleaned)
    rownames(new_cases) <- as.character(data$dates[start_:j])
    colnames(new_cases) <- as.character(data$dates[start_:j])

    #CLEARN NEW CASEES
    if(dim(new_cases)[1]>days_run){
        for(i in 1:dim(new_cases)[1]){
            if(i + days_run <= dim(new_cases)[1] ){
                index=(i+days_run):dim(new_cases)[1]
                new_cases[i, index]=NA
            }
        }
    }
    result <- model(new_cases, model_parameters, prior_list)
    result$posteriror_list$n_obs <- prior_list$n_obs + days_run

    Npost <- cbind(data$dates[j],result$Npost, N_est_true[as.character(data$dates)%in%as.character(result$Npost$dates)])
    colnames(Npost)[c(1,7)] <- c("State",'Truth')
    saveRDS(Npost, file = paste(store_data_folder,"/Npost_",j,'.rds',sep=""))
    saveRDS(result$posteriror_list, file = paste(store_data_folder,"/prior_",j,'.rds',sep=""))
}
parallel::stopCluster(cl)


