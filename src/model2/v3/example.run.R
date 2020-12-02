
graphics.off()
source(file.path("src","model2","v2","prob_dens.R"))
source(file.path("src","model2","v2","regression.R"))
source(file.path("src","model2","v2","model.R"))


# todo addapt sample N
##
sim <- 40000
data <- readRDS(file.path("data", "processed", "processed_data.rds"))
#j <- 32
j <- dim(data$detected)[1]
start_ = j - 31 # run the last 31 days

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
                   mu_lambda      = c(1,5),
                   Sigma_lambda   = 5*diag(2))

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
par(mfrow=c(2,2))
plot(result$posteriror_sample$sigma[,1])
hist(result$posteriror_sample$sigma[,1])
plot(result$posteriror_sample$sigma[,2])
hist(result$posteriror_sample$sigma[,2])
x11()
par(mfrow=c(2,2))
plot(result$posteriror_sample$lambda[,1])
hist(result$posteriror_sample$lambda[,1])
plot(result$posteriror_sample$lambda[,2])
hist(result$posteriror_sample$lambda[,2])


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

        for(j in 1:length(ind_))
                I[i,ind_[j]] <- j

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
