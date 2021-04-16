###
# running the crossvaldiation
#
###



source(file.path("src","util","util_uk.r"))
result <- readRDS(file.path("data", "processed", "processed_data_uk.rds"))

zero.report <- result$dates%in% as.Date(c("2021-01-26","2021-01-28","2021-03-01") )#rowSums(abs(diff(t(result$detected))),  na.rm=T)
max.days.to.report <- 30
result$report[,zero.report] = 0



#N.est <- sample.uk.deaths(result, max.days.to.report, samples = 1000)

pred <- uk.prediction(result = result,
              max.days.to.report = max.days.to.report,
              report.dates = c(as.Date("2021-02-28")))
pred.smooth <- gp.smooth(pred,
                         max.days.to.report = max.days.to.report)
plot(pred$date,pred$upp, type='l',col='blue')
lines(pred$date,pred$low, type='l',col='red')

lines(pred.smooth$date,pred.smooth$upp, type='l',col='green')
lines(pred.smooth$date,pred.smooth$low, type='l',col='green')
