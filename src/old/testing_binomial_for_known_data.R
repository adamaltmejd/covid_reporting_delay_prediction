unique.days = 4
maxusage.day=20
result <- readRDS(file.path("data", "processed", "processed_data.rds"))
Reported_T = result$detected
N_T <- dim(Reported_T)[1]
fixed = N_T-20
deaths_est_T <- apply(Reported_T, 1, max, na.rm=T)
deaths_est_T <- deaths_est_T[1:fixed]
data_T <- newDeaths(deaths_est_T,
                    Reported_T[1:fixed,1:fixed],
                    maxusage.day =maxusage.day)
death.remain = data_T$death.remain
report.new   = data_T$report.new
X_T <- setup_data(fixed, maxusage.day, result$dates_report[1:fixed], unique.days)

index <- upper.tri(death.remain,diag = T)
y = report.new[index]
n = death.remain[index]
index = is.na(y)==F
X <- as.matrix(X_T[index,])
y <- y[index]
n <- n[index]

fitShort <- glm( cbind(y, n) ~ X,
                 family = "binomial")
print(summary(fitShort))

