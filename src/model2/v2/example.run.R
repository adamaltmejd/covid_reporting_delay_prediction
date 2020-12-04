
graphics.off()
source(file.path("src","model2","v2","prob_dens.R"))
source(file.path("src","model2","v2","regression.R"))
source(file.path("src","model2","v2","model.R"))


days_run <- 30
# todo addapt sample N
##
sim <- 10000
data <- readRDS(file.path("data", "processed", "processed_data.rds"))
#j <- 32
j <- dim(data$detected)[1]
j <- days_run+14
start_ = j - days_run # run the last 31 days

report_cleaned <- report_clean(data$detected[start_:j,start_:j],data$dates[start_:j])
new_cases <- newCases(report_cleaned)
rownames(new_cases) <- as.character(data$dates[start_:j])
colnames(new_cases) <- as.character(data$dates[start_:j])

model_parameters <- list(sim           = sim,
                         burnin        = ceiling(0.5*sim),
                         N.days.fixed  =  3,
                         quantile      = c(0.05,0.95))

prior_list <- list(mu_beta        = c(0,0,0),
                   Sigma_beta     = 1/2*diag(3),
                   a_sigma        = c(3,3),
                   b_sigma        = c(5/2,5/2),
                   mu_lambda      = 5,
                   Sigma_lambda   = 5)

result <- model(new_cases, model_parameters, prior_list)

if(0){
    x11()
    par(mfrow=c(3,2))
    plot(result$posteriror_sample$Beta[,1])
    hist(result$posteriror_sample$Beta[,1])
    plot(result$posteriror_sample$Beta[,2])
    hist(result$posteriror_sample$Beta[,2])
    plot(result$posteriror_sample$Beta[,3])
    hist(result$posteriror_sample$Beta[,3])

    x11()
    par(mfrow=c(2,2))
    plot(result$posteriror_sample$sigma[,1])
    hist(result$posteriror_sample$sigma[,1])
    plot(result$posteriror_sample$sigma[,2])
    hist(result$posteriror_sample$sigma[,2])
}
x11()
par(mfrow=c(1,2))
plot(result$posteriror_sample$lambda)
hist(result$posteriror_sample$lambda)

if(0){
##
# analyzing output times
##
k <- sqrt(length(result$posteriror_sample$ProbMatrix[1,]))

PM <- matrix(result$posteriror_sample$ProbMatrix[1,],k,k)
I = matrix(0, nrow=k,ncol=k)
for(i in 1:k){
    index = is.na(PM[i,])==F
    if(sum(index)>0){
        ind_ <- sort(which(index))

        for(l in 1:length(ind_))
                I[i,ind_[l]] <- l

        if(ind_[1] > i+1)
            I[i, ind_[1]] = 1.5
    }
}

x11()
par(mfrow=c(3,2))
plot(colMeans(result$posteriror_sample$ProbMatrix[,c(I==1.5)]), main='t_0')
plot(colMeans(result$posteriror_sample$ProbMatrix[,c(I==1)]),  main='t_0 , no holiday')
plot(colMeans(result$posteriror_sample$ProbMatrix[,c(I==2)]), main='t_1')
plot(colMeans(result$posteriror_sample$ProbMatrix[,c(I==3)]), main='t_2')
plot(colMeans(result$posteriror_sample$ProbMatrix[,c(I==4)]), main='t_3')
cat('sd1 = ',sd(colMeans(result$posteriror_sample$ProbMatrix[,c(I==1)])),
    'm1  = ',mean(colMeans(result$posteriror_sample$ProbMatrix[,c(I==1)])),
    '\n')
cat('sd1.5 = ',sd(colMeans(result$posteriror_sample$ProbMatrix[,c(I==1.5)])),
    'm1.5  = ',mean(colMeans(result$posteriror_sample$ProbMatrix[,c(I==1.5)])),
    '\n')
cat('sd2 = ',sd(colMeans(result$posteriror_sample$ProbMatrix[,c(I==2)])),
    'm2  = ',mean(colMeans(result$posteriror_sample$ProbMatrix[,c(I==2)])),
    '\n')
cat('sd3 = ',sd(colMeans(result$posteriror_sample$ProbMatrix[,c(I==3)])),
    'm3  = ',mean(colMeans(result$posteriror_sample$ProbMatrix[,c(I==3)])),
    '\n')
cat('sd4 = ',sd(colMeans(result$posteriror_sample$ProbMatrix[,c(I==4)])),
    'm4  = ',mean(colMeans(result$posteriror_sample$ProbMatrix[,c(I==4)])),
    '\n')

cat('sd5 = ',sd(colMeans(result$posteriror_sample$ProbMatrix[,c(I==5)])),
    'm5  = ',mean(colMeans(result$posteriror_sample$ProbMatrix[,c(I==5)])),
    '\n')
cat('sd6 = ',sd(colMeans(result$posteriror_sample$ProbMatrix[,c(I==6)])),
    'm6  = ',mean(colMeans(result$posteriror_sample$ProbMatrix[,c(I==6)])),
    '\n')
cat('sd7 = ',sd(colMeans(result$posteriror_sample$ProbMatrix[,c(I==7)])),
    'm7  = ',mean(colMeans(result$posteriror_sample$ProbMatrix[,c(I==7)])),
    '\n')
}
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
