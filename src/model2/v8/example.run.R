
graphics.off()
source(file.path("src","model2","v8","prob_dens.R"))
source(file.path("src","model2","v8","regression.R"))
source(file.path("src","model2","v8","model.R"))


days_run <- 30
# todo addapt sample N
##
sim <- 10000
data <- readRDS(file.path("data", "processed", "processed_data.rds"))
j <-  55
#j <- dim(data$detected)[1]-10
start_ = j - days_run # run the last 31 days

report_cleaned <- report_clean(data$detected[start_:j,start_:j],data$dates[start_:j])
new_cases <- newCases(report_cleaned)
rownames(new_cases) <- as.character(data$dates[start_:j])
colnames(new_cases) <- as.character(data$dates[start_:j])

model_parameters <- list(sim           = sim,
                         burnin        = ceiling(0.5*sim),
                         N.days.fixed  =  3,
                         quantile      = c(0.05,0.95))

prior_list <- list(mu_beta        = c(0,0,0,0),
                   Sigma_beta     = 2*diag(4),
                   a_sigma        = c(3,3),
                   b_sigma        = c(5/2,5/2),
                   mu_GP          = c(0,10),
                   mu_theta_GP    = c(0,-2),
                   Sigma_theta_GP = diag(c(1,1)),
                   mu_phi         = 1,
                   sigma_phi      = 1)

result <- model(new_cases, model_parameters, prior_list)

x11()
par(mfrow=c(4,2))
plot(result$posteriror_sample$Beta[,1])
hist(result$posteriror_sample$Beta[,1])
plot(result$posteriror_sample$Beta[,2])
hist(result$posteriror_sample$Beta[,2])
plot(result$posteriror_sample$Beta[,3])
hist(result$posteriror_sample$Beta[,3])
plot(result$posteriror_sample$Beta[,4])
hist(result$posteriror_sample$Beta[,4])



x11()
par(mfrow=c(2,1))
m<-length(as.Date(result$Npost$dates))
ylim <- c(min(c(colMeans(exp(result$posteriror_sample$theta)),result$Npost$lCI)),
          max(c(colMeans(exp(result$posteriror_sample$theta)),result$Npost$uCI)))
plot(as.Date(result$Npost$dates),colMeans(exp(result$posteriror_sample$theta)), type='l',ylim=ylim)
points(as.Date(result$Npost$dates),apply(data$detected,1, max, na.rm=T)[start_:(m+start_-1)], col='red')
lines(as.Date(result$Npost$dates),result$Npost$median,col='blue')
lines(as.Date(result$Npost$dates),result$Npost$lCI,col='blue')
lines(as.Date(result$Npost$dates),result$Npost$uCI,col='blue')
plot(result$posteriror_sample$theta[,dim(result$posteriror_sample$theta)[2]])
print(cbind(data$dates[j]-as.Date(result$Npost$dates),round(result$Npost$uCI-result$Npost$lCI,2)))
x11()
par(mfrow=c(2,1))
plot(result$posteriror_sample$phi)
hist(result$posteriror_sample$phi, prob=T)
curve(dlnorm(x, meanlog=result$posteriror_list$mu_phi, sdlog=result$posteriror_list$sigma_phi),
      col="darkblue", lwd=2, add=TRUE, yaxt="n")
x11()
par(mfrow=c(2,2))
plot(result$posteriror_sample$sigma_theta^2)
hist(result$posteriror_sample$sigma_theta^2, prob=T)
curve(dlnorm(x, meanlog=result$posteriror_list$mu_theta_GP[2], sdlog=sqrt(result$posteriror_list$Sigma_theta_GP[2,2])),
      col="darkblue", lwd=2, add=TRUE, yaxt="n")
plot(result$posteriror_sample$kappa)
hist(result$posteriror_sample$kappa, prob=T)
curve(dlnorm(x, meanlog=result$posteriror_list$mu_theta_GP[1], sdlog=sqrt(result$posteriror_list$Sigma_theta_GP[1,1])),
      col="darkblue", lwd=2, add=TRUE, yaxt="n")

x11()
par(mfrow=c(1,1))
plot(as.Date(result$Npost$dates),result$Npost$median,col='blue',type='l')
end_ <- start_ + length(result$Npost$median) -1
points(as.Date(result$Npost$dates),apply(data$detected,1, max, na.rm=T)[start_:end_], col='red')

lines(as.Date(result$Npost$dates),result$Npost$lCI,col='blue')
lines(as.Date(result$Npost$dates),result$Npost$uCI,col='blue')
print(cbind(data$dates[j]-as.Date(result$Npost$dates),
            result$Npost$median,apply(data$detected,1, max, na.rm=T)[start_:end_],
            round(result$Npost$uCI-result$Npost$lCI,2)))


