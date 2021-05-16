###
# running the crossvaldiation
#
###



source(file.path("src","util","util_swe.r"))
result <- readRDS(file.path("data", "processed", "processed_data_sweden.rds"))


#remove data before  2020-07-01
index <- result$dates_report >= as.Date("2020-07-01")
result$detected             = result$detected[index,index]
result$report               = result$report[index,index]
result$dates                = result$dates[index]
result$dates_report         = result$dates_report[index]
result$dates_not_reported   = result$dates_not_reported[index]
max.days.to.report <- 30


target <- data.frame(reported = result$detected[row(result$detected)+max.days.to.report==col(result$detected)])
target$dates <- result$dates_report[1:length(target$reported)]
#remove all index above max
result$detected[row(result$detected)+max.days.to.report<col(result$detected)]=NA
#N.est <- sample.uk.deaths(result, max.days.to.report, samples = 1000)

pred <- swe.prediction(result = result,
              max.days.to.report = max.days.to.report,
              report.dates = c(as.Date("2021-03-20")),
              target= target)
pred.smooth <- gp.smooth(pred,
                         max.days.to.report = max.days.to.report)
plot(pred$state,pred$ci_lower, type='l',col='blue',ylim=c(0,140))
lines(pred$state,pred$ci_upper, type='l',col='red')

lines(pred.smooth$state,pred.smooth$ci_lower, type='l',col='green')
lines(pred.smooth$state,pred.smooth$ci_upper, type='l',col='green')
points(target$dates,target$reported)
