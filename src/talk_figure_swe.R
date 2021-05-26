save.fig=T
smooth = T
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
points(swe_data$dates[(i-30):i],swe_data$detected[(i-30):i,i],col='black',pch=20)

if(save.fig)
    dev.off()

if(save.fig)
    pdf("swe_crps.pdf")

plot(dts_smooth[state==date-4,state],dts_smooth[state==date-4,CRPS],xlab='date',ylab='CRPS')

if(save.fig)
    dev.off()

dts_uk <-read_fst( file.path("data", "model_predictions_full_UK.fst"))
dts_uk <-as.data.table(dts_uk)
dts_smooth_uk <-read_fst( file.path("data", "model_predictions_full_UK_smooth.fst"))
dts_smooth_uk <-as.data.table(dts_smooth_uk)


date_plots <- as.Date("2021-01-21")+0:10

if(smooth){
    dts_plot <- dts_smooth_uk
    file.name <- "uk_smooth_"
}else{
    dts_plot <- dts_uk
    file.name <- "uk_"
}
for(j in 1:length(date_plots)){
    date_plot = date_plots[j]
    if(save.fig)
        pdf(paste(file.name,j,'.pdf',sep=''))
    plot(dts_plot[date==date_plot,state],dts_plot[date==date_plot,target],ylim=c(0,1600), xlab='date',ylab='dea')
    lines(dts_plot[date==date_plot,state],dts_plot[date==date_plot,ci_upper],col='blue')
    lines(dts_plot[date==date_plot,state],dts_plot[date==date_plot,ci_lower],col='red')
    uk_data <- readRDS(file.path("data", "processed", "processed_data_uk.rds"))
    i = which(uk_data$dates==date_plot)
    points(uk_data$dates[(i-30):i],uk_data$detected[(i-30):i,i],col='black',pch=20)

    if(save.fig)
        dev.off()
}
