###
#  Building the data file for the model
#  main thing done is averaing over overcount of deaths and creates and average of them see line:
# 36-48
# load files:
# data/processed/deaths_dt.fst (data_processing.R)
#
# generates the files:
# data/processed/processed_data.rds
#
# detected is an NxN matrix
# detected[i,]  - represtents how many reported dead for dates[i],
#                 the columns represents reports days
# detected[,i]  - reprsentes reported deaths at rerport dayes dates[i]
##
library(Matrix)
library(stringr)
library(dplyr)
library(fst)


buildData <- function(country = "sweden"){

    if(country == "sweden"){
        deaths_dt <- read_fst(file.path("data", "processed", "deaths_dt.fst"), as.data.table = TRUE)
        deaths_dt <- deaths_dt[date > "2020-04-01", .(date, publication_date, N)]
    }else if(country=="uk"){
        deaths_dt <- read_fst(file.path("data", "processed", "deaths_dt_UK.fst"), as.data.table = T)
        deaths_dt <- deaths_dt[, .(date, publication_date, N)]
        deaths_dt <- deaths_dt[date >= "2020-09-01", .(date, publication_date, N)]

    }
    ##
    #only relevant stuff below
    ##
    res2 <- deaths_dt %>% tidyr::spread(publication_date, N)

    # dayes with no reporting
    index_NoReport <- res2[,1]$date%in%unique(deaths_dt$publication_date)==F


    # create a detected matrix for date x date
    n.days <- dim(res2)[1]
    detected <- matrix(NA, nrow= n.days,
                          ncol =n.days)
    detected[,index_NoReport == F] = as.matrix(res2[,2:dim(res2)[2]])
    colnames(detected) <- as.character((res2[,1]$date))
    rownames(detected) <- as.character((res2[,1]$date))


    repeated = 1
    while(repeated != 0){
        repeated = 0
        for(i in 1:(dim(detected)[1]-1)){
            for(j in i:(dim(detected)[2]-1)){
                if(is.na(detected[i,j])==F & is.na(detected[i,j+1])==F){
                    if(detected[i,j] > detected[i,j+1]){
                        temp = ceiling(0.5*detected[i,j] + 0.5*detected[i,j+1])
                        detected[i,j] =temp
                        detected[i,j+1] =temp
                        repeated = 1
                    }
                }
            }
        }
    }
    result <- list(detected = detected,
                   dates = res2[,1]$date,
                   dates_report=unique(deaths_dt$publication_date),
                   dates_not_reported = index_NoReport)
    if(country == "sweden"){
        saveRDS(result, file.path("data", "processed", "processed_data.rds"))
    }else if(country=="uk"){
        saveRDS(result, file.path("data", "processed", "processed_data_uk.rds"))
    }
}

buildData("sweden")
buildData("uk")
