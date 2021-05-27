###
# running the crossvaldiation
#
###



source(file.path("src","util","util_uk.r"))
result <- readRDS(file.path("data", "processed", "processed_data_uk.rds"))

#zero.report <- result$dates%in% as.Date(c("2021-01-26","2021-01-28","2021-03-01") )#rowSums(abs(diff(t(result$detected))),  na.rm=T)
max.days.to.report <- 30
#result$report[,zero.report] = 0

target <- data.frame(reported = result$detected[row(result$detected)+max.days.to.report==col(result$detected)])
target$dates <- result$dates_report[1:length(target$reported)]
#remove all index above max
result$detected[row(result$detected)+max.days.to.report<col(result$detected)]=NA
#N.est <- sample.uk.deaths(result, max.days.to.report, samples = 1000)

pred <- uk.prediction(result = result,
              max.days.to.report = max.days.to.report,
              report.dates = c(as.Date("2021-01-12")),
              target= target)
pred.smooth <- gp.smooth(pred,
                         max.days.to.report = max.days.to.report)
plot(pred$date,pred$target,ylim=c(0,max(pred.smooth$ci_upper)))
lines(pred$date,pred$predicted_deaths, type='l',col='blue')
lines(pred$date,pred$ci_upper, type='l',col='red')
lines(pred$date,pred$ci_lower, type='l',col='red')

lines(pred.smooth$date,pred.smooth$ci_upper, type='l',col='green')
lines(pred.smooth$date,pred.smooth$ci_lower, type='l',col='green')
#2021-01-25
