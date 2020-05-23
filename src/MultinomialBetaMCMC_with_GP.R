##
# fitting maximum likelhiood multinomial using logit p
# then profile likelihood of N
##
library(data.table)
source(file.path("src", "util.r"))
source(file.path("src", "MH.R"))
source(file.path("src", "functions.R"))
source(file.path("src", "GPutil.R"))
source(file.path("src", "modelbenchmarkutil.R"))
set.seed(4)

path.to.files <- file.path("..","data")
save.file = F
MCMC_sim <- 10000
burnin_p = 0.33
deaths_sim <- 5
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

day <- 16
res_save<-benchmark_BetaGP_j(day,
                             result,true.day,
                             maxusage.day,
                             unique.days,
                             MCMC_sim,
                             burnin_p,
                             deaths_sim = deaths_sim,
                             prior = c(0,0))

CI <-apply(res_save$Death_est,2 , function(x){ quantile(x,c(0.05,0.5,0.95))})

out <- data.table(date = result$dates[1:day],
                      ci_upper  = CI[3,],
                      ci_lower  = CI[1,],
                      predicted_deaths = CI[2,])
plot(out$date, out$predicted_deaths, type='l', ylim = c(0,max(out$ci_upper)),col='blue')
points(out$date, result$detected[1:day,day])
points(out$date, deaths_est[1:day], pch=3)
lines(out$date, out$ci_upper,  ylim = c(0,max(out$ci_upper)),col='red')
lines(out$date, out$ci_lower, ylim = c(0,max(out$ci_upper)),col='red')
if(save.file)
    write_fst(out, file.path("data", "processed", "model_predict.fst"))


