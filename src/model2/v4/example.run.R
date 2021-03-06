
graphics.off()
source(file.path("src","model2","v4","prob_dens.R"))
source(file.path("src","model2","v4","regression.R"))
source(file.path("src","model2","v4","model.R"))


days_run <- 30
# todo addapt sample N
##
sim <- 1000
data <- readRDS(file.path("data", "processed", "processed_data.rds"))
j <- days_run+14
#j <- dim(data$detected)[1]-10
start_ = j - days_run # run the last 31 days

report_cleaned <- report_clean(data$detected[start_:j,start_:j],data$dates[start_:j])
new_cases <- newCases(report_cleaned)
rownames(new_cases) <- as.character(data$dates[start_:j])
colnames(new_cases) <- as.character(data$dates[start_:j])

model_parameters <- list(sim           = sim,
                         burnin        = ceiling(0.5*sim),
                         N.days.fixed  =  3,
                         quantile      = c(0.1,0.9))

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

result <- model(new_cases, model_parameters, prior_list)

x11()
par(mfrow=c(3,2))
plot(result$posteriror_sample$Beta[,1])
hist(result$posteriror_sample$Beta[,1])
plot(result$posteriror_sample$Beta[,2])
hist(result$posteriror_sample$Beta[,2])
plot(result$posteriror_sample$Beta[,3])
hist(result$posteriror_sample$Beta[,3])


x11()
par(mfrow=c(2,1))
plot(result$posteriror_sample$lambda)
hist(result$posteriror_sample$lambda)

x11()
par(mfrow=c(2,1))
ylim <- c(min(c(colMeans(exp(result$posteriror_sample$theta)),result$Npost$lCI)),
          max(c(colMeans(exp(result$posteriror_sample$theta)),result$Npost$uCI)))
plot(colMeans(exp(result$posteriror_sample$theta)), type='l',ylim=ylim)
points(apply(data$detected,1, max, na.rm=T)[start_:j], col='red')
lines(result$Npost$median,col='blue')
lines(result$Npost$lCI,col='blue')
lines(result$Npost$uCI,col='blue')
plot(result$posteriror_sample$theta[,dim(result$posteriror_sample$theta)[2]])
print(cbind(data$dates[j]-as.Date(result$Npost$dates),round(result$Npost$uCI-result$Npost$lCI,2)))
x11()
par(mfrow=c(2,1))
plot(result$posteriror_sample$phi)
hist(result$posteriror_sample$phi, prob=T)
curve(dlnorm(x, meanlog=result$posteriror_list$mu_phi, sdlog=result$posteriror_list$sigma_phi),
      col="darkblue", lwd=2, add=TRUE, yaxt="n")
x11()
par(mfrow=c(2,1))
plot(result$posteriror_sample$sigma_theta)
hist(result$posteriror_sample$sigma_theta^2, prob=T)
curve(dinvgamma(x, shape=result$posteriror_list$a_sigma_theta, rate=result$posteriror_list$b_sigma_theta),
      col="darkblue", lwd=2, add=TRUE, yaxt="n")

x11()
par(mfrow=c(1,1))
plot(result$Npost$median,col='blue',type='l')
end_ <- start_ + length(result$Npost$median) -1
points(apply(data$detected,1, max, na.rm=T)[start_:end_], col='red')

lines(result$Npost$lCI,col='blue')
lines(result$Npost$uCI,col='blue')
print(cbind(data$dates[j]-as.Date(result$Npost$dates),
            result$Npost$median,apply(data$detected,1, max, na.rm=T)[start_:end_],
            round(result$Npost$uCI-result$Npost$lCI,2)))


