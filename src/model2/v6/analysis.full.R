##
# first run full.run to get files
#
#
##
graphics.off()
store_data_folder <- file.path("data","tmp","model2","v6")
data_files <- list.files(store_data_folder,pattern='Npost*')

data <- c()
for(file in data_files){
    data <- rbind(data, readRDS(paste(store_data_folder,'/',file,sep="")))
}
data$df <- as.Date(data$State)-as.Date(data$dates)

days = 1:14
for(day in days){
    index <- data$df==day
    data_temp <- data[index,]
    CI_cov <- (data_temp$Truth >= data_temp$lCI) &(data_temp$Truth <= data_temp$uCI)

    cat(' day =  ',day,': ')
    cat(' mabs = ', round(mean(abs(data_temp$Truth-data_temp$median)),2))
    cat(' rmse = ',round(sqrt(mean((data_temp$Truth-data_temp$median)^2)),2))
    cat(' CI =  ',round(mean(CI_cov),2))
    cat('\n')
}
index <- data$df==7
data_temp <- data[index,]
CI_cov <- (data_temp$Truth >= data_temp$lCI) &(data_temp$Truth <= data_temp$uCI)

x11()
par(mfrow=c(2,2))
plot(1:length(CI_cov),cumsum(CI_cov)/(1:length(CI_cov)),xlab='days',ylab='coverage',main='7 days')
plot(data_temp$Truth,CI_cov)
plot(data_temp$Truth,data_temp$Truth-data_temp$median)
