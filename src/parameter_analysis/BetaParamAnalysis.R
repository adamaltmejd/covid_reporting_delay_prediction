##
#  used for analazing the convergences of the probability of detection parameters
#  i.e. beta distribution
#  load file:
#  data/processed/processed_data.rds (buildData.R)

#  generates the files:
##
library(tidyr)
source(file.path("src", "util","util.r"))
source(file.path("src", "util","MH.R"))
source(file.path("src", "util","functions.R"))
source(file.path("src", "util","GPutil.R"))
source(file.path("src","util","MLbeta.R"))
lag <- 1
days_to_pred <- 40

predicition.list <-list()
result <- readRDS(file.path("data", "processed", "processed_data.rds"))
npar = 2



j <- 30
param_alpha <- c()
param_beta  <- c()
for(j in  31:length(result$dates)){
start_  <- max(1,j-days_to_pred)
result_j <- result

result_j$detected     <- result_j$detected[start_:j,start_:j]
result_j$dates        <- result_j$dates[start_:j]
result_j$dates_report <-  result_j$dates_report[start_:j]
result_j$dates_not_reported <- result_j$dates_not_reported[start_:j]
res <-  ML_betaBin(result_j,
                   lag,
                   npar)
param_alpha <-cbind(param_alpha, res$alpha_X)
param_beta <-cbind(param_beta, res$beta_X)
}
