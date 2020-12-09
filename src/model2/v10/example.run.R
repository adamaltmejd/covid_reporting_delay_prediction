
graphics.off()
source(file.path("src","model2","v10","prob_dens.R"))
source(file.path("src","model2","v10","regression.R"))
source(file.path("src","model2","v10","model.R"))


days_run <- 30
# todo addapt sample N
##
sim <- 40000
data <- readRDS(file.path("data", "processed", "processed_data.rds"))
j <- dim(data$detected)[1] -29
start_ = j - days_run # run the last 31 days

report_cleaned <- report_clean(data$detected[start_:j,start_:j],data$dates[start_:j])
new_cases <- newCases(report_cleaned)
rownames(new_cases) <- as.character(data$dates[start_:j])
colnames(new_cases) <- as.character(data$dates[start_:j])

model_parameters <- list(sim           = sim,
                         burnin        = ceiling(0.5*sim),
                         N.days.fixed  =  3,
                         quantile      = c(0.025,0.975))

prior_list <- list(mu_beta        = c(0,0,0,0,0),
                   Sigma_beta     = 50*diag(5),
                   a_sigma        = c(3,3,3),
                   b_sigma        = c(10/2,10/2,10/2),
                   theta_mu         = c(0,0),
                   theta_Sigma    =  5*diag(2),
                   mu_phi         = 1,
                   sigma_phi      = 1)

result <- model(new_cases, model_parameters, prior_list)

x11()
par(mfrow=c(3,2))
plot(result$posteriror_sample$Beta[,1],main='first')
hist(result$posteriror_sample$Beta[,1])
plot(result$posteriror_sample$Beta[,2],main='second,third, fourth')
hist(result$posteriror_sample$Beta[,2])
plot(result$posteriror_sample$Beta[,3],main='rest')
hist(result$posteriror_sample$Beta[,3])
x11()
par(mfrow=c(2,2))
plot(result$posteriror_sample$Beta[,4],main='Tuesday')
hist(result$posteriror_sample$Beta[,4])
plot(result$posteriror_sample$Beta[,5],main='Wednseday')
hist(result$posteriror_sample$Beta[,5])




x11()
par(mfrow=c(2,1))
plot(result$posteriror_sample$phi)
hist(result$posteriror_sample$phi, prob=T)
curve(dlnorm(x, meanlog=result$posteriror_list$mu_phi, sdlog=result$posteriror_list$sigma_phi),
      col="darkblue", lwd=2, add=TRUE, yaxt="n")
x11()
par(mfrow=c(2,1))
plot(result$posteriror_sample$theta[,1])
plot(result$posteriror_sample$theta[,2])
x11()
par(mfrow=c(1,1))
m<-length(as.Date(result$Npost$dates))
truth <- apply(data$detected,1, max, na.rm=T)[start_:(m+start_-1)]

ylim <- c(min(c(colMeans(exp(result$posteriror_sample$theta)),result$Npost$lCI, truth)),
          max(c(colMeans(exp(result$posteriror_sample$theta)),result$Npost$uCI, truth)))
plot(as.Date(result$Npost$dates),
     truth,
     col='red',
     ylim=ylim,
     type='p')
lines(as.Date(result$Npost$dates),result$Npost$median,col='blue')
lines(as.Date(result$Npost$dates),result$Npost$lCI,col='blue')
lines(as.Date(result$Npost$dates),result$Npost$uCI,col='blue')

Res <- colMeans(result$posteriror_sample$t)
Prob_mat <- matrix(result$posteriror_sample$ProbMatrix[1,], nrow= sqrt(length(result$posteriror_sample$ProbMatrix[1,])))
print(Prob_mat)
