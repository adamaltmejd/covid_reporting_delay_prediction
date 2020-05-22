##
# fitting maximum likelhiood multinomial using logit p
# then profile likelihood of N
##
set.seed(2)
library(data.table)
source(file.path("src", "util.r"))
source(file.path("src", "MH.R"))
source(file.path("src", "functions.R"))
source(file.path("src", "GPutil.R"))
source(file.path("src", "modelbenchmarkutil.R"))


path.to.files <- file.path("..","data")
MCMC_sim <- 10000
burnin_p = 0.33
deaths_sim <- 10
maxusage.day = 20 #must be less then N
unique.days  = 7
true.day = 5
result <- readRDS(file.path("data", "processed", "processed_data.rds"))

Reported = result$detected
N <- dim(Reported)[1]

deaths_est <- apply(Reported, 1, max, na.rm=T)

data_ <- newDeaths(deaths_est,
                   Reported,
                   maxusage.day =maxusage.day)

res_save<-benchmark_BetaGP_j(N,
                             result,true.day,
                             maxusage.day,
                             unique.days,
                             MCMC_sim,
                             burnin_p,
                             deaths_sim = deaths_sim,
                             prior = c(0,0))

CI <-apply(res_save$Death_est,2 , function(x){ quantile(x,c(0.05,0.5,0.95))})

out <- data.table(date = result$dates,
                      ci_upper  = CI[1,],
                      ci_lower  = CI[3,],
                      predicted_deaths = CI[2,])
write_fst(out, file.path("data", "processed", "model_predict.fst"))
