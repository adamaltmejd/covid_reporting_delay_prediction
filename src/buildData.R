###
# load data
##
library(Matrix)
library(stringr)
library(dplyr)
library(fst)

# source("src/functions.R")

# Sys.setlocale("LC_TIME", 'en_US.UTF-8')
# #download_latest_fhm()
# path.to.files <- file.path("data")
# fhm_files = list_fhm_files(folder = paste(path.to.files,"/FHM/",sep=""))
# death_dts <- c()
# for(i in 1:length(fhm_files)){
#   date_in <- as.Date(fhm_files[i]%>% str_match_all("[0-9]{4}-[0-9]{2}-[0-9]{2}") %>%unlist)
#   if(date_in > as.Date('2020-04-01'))
#     death_dts = rbind(death_dts,load_fhm(fhm_files[i]))
# }
#
# # res <- death_dts %>% dplyr::group_by(date, publication_date) %>% dplyr::spread(publication_date, value=N)
# res <- death_dts %>% dplyr::group_by(date, publication_date)
# death_dt <- data.table(death_dts)
# setkey(death_dt, publication_date, date)
# death_dt <-death_dt[!is.na(date) & publication_date > "2020-04-02" & date > "2020-04-02"]

deaths_dt <- read_fst(file.path("data", "processed", "deaths_dt.fst"), as.data.table = TRUE)
deaths_dt <- deaths_dt[date > "2020-04-01", .(date, publication_date, N)]
##
#only relevant stuff below
##
res2 <- deaths_dt %>% tidyr::spread(publication_date, N)

deaths_dt[, .N, date]
deaths_dt[is.na(date), .N]
deaths_dt[date < "2020-03-11"]

detected <- as.matrix(res2[,2:dim(res2)[2]])
detected <- cbind(matrix(NA,dim(detected)[1],dim(detected)[1]-dim(detected)[2] ) , detected)
detected[lower.tri(detected)]=NA
colnames(detected) <- c(res2[,1]$date)
d_ <- diag(detected)
d_[is.na(d_)] <- 0
diag(detected) <- d_

result <- list(detected = detected, dates = res2[,1]$date, dates_report=unique(deaths_dt$publication_date))
saveRDS(result, file.path("data", "processed", "processed_data.rds"))

