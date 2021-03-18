
graphics.off()
source(file.path("src","model2","v5","prob_dens.R"))
source(file.path("src","model2","v5","regression.R"))
source(file.path("src","model2","v5","model.R"))

graphics.off()
data <- readRDS(file.path("data", "processed", "processed_data.rds"))
#data <- readRDS(file.path("data", "processed", "processed_data_icu.rds"))


report_cleaned <- report_clean(data$detected,data$dates)
new_cases <- newCases(report_cleaned)
rownames(new_cases) <- as.character(data$dates)
colnames(new_cases) <- as.character(data$dates)
diag(new_cases) <- NA
N <- apply(new_cases,1,sum, na.rm=T)
lag <- 15
to. <- 30
Prob <- matrix(0,nrow=dim(new_cases)[1]-to.,ncol=lag)
for(i in 1:(dim(new_cases)[1]-to.)){
  cases <- new_cases[i,]
  m <- min(length(cases)-i,lag)
  cases <- cases[i+ (1:m)]
  cases <- cases[is.na(cases)==F ]

  Prob[i,1:m] = cases[1:m]/sum(new_cases[i,], na.rm=T)
  #Prob[i,1:m] = cases[1:m]/N[i]
}
dates <- as.Date(data$dates[1:(dim(new_cases)[1]-to.)])
N <- N[1:(dim(new_cases)[1]-to.)]
N_min <- 6
x11()
par(mfrow=c(2,2))
plot(dates[N>N_min],apply(Prob[N>N_min,1:3,drop=F],1,sum))
plot(dates[N>N_min],apply(Prob[N>N_min,1:5,drop=F],1,sum))
plot(dates[N>N_min],apply(Prob[N>N_min,1:6,drop=F],1,sum))
plot(dates[N>N_min],apply(Prob[N>N_min,1:7, drop=F],1,sum))
if(0){
x11()
index <- weekdays(data$dates)%in%"Tuesday"
D     <- apply(new_cases,2,sum,na.rm=T)
plot(data$dates[index],D[index],col='red',type='l')
index <- weekdays(data$dates)%in%"Wednesday"
lines(data$dates[index],D[index],col='blue')
index <- weekdays(data$dates)%in%"Thursday"
lines(data$dates[index],D[index],col='green')
index <- weekdays(data$dates)%in%"Friday"
lines(data$dates[index],D[index],col='black')
}


