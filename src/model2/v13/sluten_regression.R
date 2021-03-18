##
# analysis of slutenvard
#
##
library(gamm4)

source(file.path("src","model2","v12","icu_regression.R"))
lag <- 7
deaths <- readRDS(file.path("data", "processed", "processed_data.rds"))
data.vard <- read.csv(file.path("data", "vard", "data-fJECo.csv"))

n = dim(data.vard)[1]
slutenvard <- as.numeric(data.vard[,3])
intensvard <-  as.numeric(data.vard[,5])
days       <- 1:length(slutenvard)
data.list.sluten <- data.frame(days= days[!is.na(slutenvard)],
                        date = as.Date(data.vard[!is.na(slutenvard),1]),
                        sluten   = slutenvard[!is.na(slutenvard)])
model.fit.sluten <- gam(sluten~  s(days), data=data.list.sluten, family = nb(),method="REML" )
data.list.iva <- data.frame(days= days[!is.na(intensvard)],
                        date = as.Date(data.vard[!is.na(intensvard),1]),
                        iva   = intensvard[!is.na(intensvard)])
model.fit.iva <- gam(iva~  s(days), data=data.list.iva, family = nb(),method="REML" )


cov.data <- data.frame( date  =as.Date(data.vard[,1]),
                        theta_sluten = predict(model.fit.sluten,newdata = data.frame(days = days)),
                        theta_iva = predict(model.fit.iva,newdata = data.frame(days = days)))
Y <- apply(deaths$detected,1,max,na.rm=T)
death.list = data.frame( date = as.Date(deaths$dates), death = Y)

lags <- 1:20
lags.iva    <- rep(NA,length(lags))
lags.sluten <- rep(NA,length(lags))
for(i in 1:length(lags)){
    lag <- lags[i]
    cov.data.lag = data.frame( date         = c(cov.data$date, cov.data$date[n]+1:lag),
                               theta_sluten = c(rep(0,lag),cov.data$theta_sluten),
                               theta_iva    = c(rep(0,lag),cov.data$theta_iva))

    data.total <- merge(cov.data.lag,death.list, by ="date")

    data.total.rem <- data.total[30:(dim(data.total)[1]-15),]

    model.fit2.iva <- gam(death~  theta_iva, data=data.total.rem, family = nb(),method="REML" )
    model.fit2.sluten <- gam(death~  theta_sluten , data=data.total.rem, family = nb(),method="REML" )
    model.fit2.joint <- gam(death~  theta_sluten  + theta_iva, data=data.total.rem, family = nb(),method="REML" )
    lags.iva[i]    <- -model.fit2.iva$aic
    lags.sluten[i] <- -model.fit2.sluten$aic
    if(-model.fit2.iva$aic == max(lags.iva, na.rm=T))
        model.iva <- model.fit2.iva
    if(-model.fit2.sluten$aic == max(lags.sluten, na.rm=T))
        model.sluten <- model.fit2.sluten


}
plot(lags, lags.iva)
points(lags, lags.sluten,col='red')
x11()
plot(data.total.rem$date,data.total.rem$death, xlab='date',ylab='Y')
lines(data.total.rem$date,exp(predict(model.iva)),col='red')
lines(data.total.rem$date,exp(predict(model.sluten)),col='blue')
legend(data.total.rem$date[110], 80, legend=c("iva pred lag = 6 days", "sluten lag= 7 days"),
       col=c("red", "blue"), lty=1:1, cex=0.5)
