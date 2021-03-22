##
# UK analysis
#
#
##

source(file.path("src","util","util.r"))
result <- readRDS(file.path("data", "processed", "processed_data_uk.rds"))

max.days.to.report <- 30

# assume cases are all reported in 30 days
N <- deaths_at_t(result$detected,max.days.to.report)

fixed.dates <- sum(is.na(N)==F)
N.fixed <- N[1:fixed.dates]
Data <- newDeaths(N.fixed, reports = result$detected[1:fixed.dates,1:fixed.dates ], maxusage.day = max.days.to.report)
# plot(diag(Data$death.remain[,-(1:8)])/N.fixed[1:(length(N.fixed)-8)])
k <- 6
plot(diag(Data$report.new[,-(1:k)])/(diag(Data$report.new[,-(1:k)])+diag(Data$death.remain[,-(1:k)])))

# 5 above everyting stable around 20% of cases reported each day


X_T <- setup_data(fixed.dates, max.days.to.report, result$dates_report[1:fixed.dates], 5)
index <- upper.tri(Data$death.remain,diag = T)
y = Data$report.new[index]
n = Data$death.remain[index]
index = is.na(y)==F
X <- as.matrix(X_T[index,])
y <- y[index]
n <- n[index]

fitShort <- glm( cbind(y, n) ~ X,
                 family = "binomial")

day2 = X[,2] == 1
y.day2 <- y[day2]>0
sunday <- X[day2,8]
sunday <- X[day2,8]
fitShort <- glm( y.day2~ 1 + sunday,
                 family = "binomial")
