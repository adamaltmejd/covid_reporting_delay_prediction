save.fig=T
smooth = F
graphics.off()
date_plot <- "2020-11-17"
dts<-read_fst( file.path("data", "model_predictions_full_SWE.fst"))
dts <-as.data.table(dts)
dts_smooth<-read_fst( file.path("data", "model_predictions_full_SWE_smooth.fst"))
dts_smooth <-as.data.table(dts_smooth)
if(smooth){
    dts_plot <- dts_smooth
    file.name <- "swe_smooth_2011117.pdf"
}else{
   dts_plot <- dts
   file.name <- "swe_2011117.pdf"
}
if(save.fig)
    pdf(file.name)
plot(dts_plot[date==date_plot,state],dts_plot[date==date_plot,target],ylim=c(0,60), xlab='date',ylab='dea')
lines(dts_plot[date==date_plot,state],dts_plot[date==date_plot,ci_upper],col='blue')
lines(dts_plot[date==date_plot,state],dts_plot[date==date_plot,ci_lower],col='red')

swe_data <- readRDS(file.path("data", "processed", "processed_data_sweden.rds"))

start_date <- as.Date("2020-07-01")
index <- swe_data$dates_report >= start_date
swe_data$detected             = swe_data$detected[index,index]
swe_data$report               = swe_data$report[index,index]
swe_data$dates                = swe_data$dates[index]
swe_data$dates_report         = swe_data$dates_report[index]
swe_data$dates_not_reported   = swe_data$dates_not_reported[index]
i = which(swe_data$dates==date_plot)
points(swe_data$dates[1:i],swe_data$detected[1:i,i],col='black',pch=20)

if(save.fig)
    dev.off()
