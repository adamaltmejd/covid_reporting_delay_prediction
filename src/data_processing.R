##
# building base data file from FHM data files
#
# reads files:
#  data/FHM/Folkhalsomyndigheten_Covid19_YYYY-MM-DD.xlsx
# (external files avialable at <https://github.com/adamaltmejd/covid>)
#
# generates the files:
# /data/processed/deaths_dt.fst
##

library(data.table)
library(readxl)
library(stringr)
library(fst)

can_be_numeric <- function(x) {
    # Check if vector can be converted to numeric
    stopifnot(is.atomic(x) || is.list(x)) # check if x is a vector
    numNAs <- sum(is.na(x))
    numNAs_new <- suppressWarnings(sum(is.na(as.numeric(x))))
    return(numNAs_new == numNAs)
}

get_record_date <- function(f) {
    sheets <- excel_sheets(f)
    ret <- as.Date(sub("^FOHM ", "", sheets[length(sheets)]), format="%d %b %Y")
    if (is.na(ret)) ret <- as.Date(sub("^FOHM ", "", sheets[grep("FOHM", sheets)]), format="%d %b %Y")
    return(ret)
}

load_fhm_deaths <- function(f) {
    require(data.table)
    require(readxl)
    require(stringr)

    date <- as.Date(str_extract(f, "[0-9]{4}-[0-9]{2}-[0-9]{2}"))
    if (date <= as.Date("2020-04-01")) return(NULL)

    DT <- data.table((
        read_excel(path = f, sheet = 2, col_types = c("text", "numeric"), progress = FALSE)
    ))

    setnames(DT, c("date", "N"))

    DT[(tolower(date) %in% c("uppgift saknas", "uppgift saknaa", "uppgift saknas+a1")), date := NA]

    if (can_be_numeric(DT[, date])) {
        DT[, date := as.Date(as.numeric(date), origin = "1899-12-30")]
    } else {
        DT[, date := as.Date(date)]
    }

    # Ensure starting point is March 1st, and that all dates have a value
    publication_date <- get_record_date(f)
    DT <- merge(DT, data.table(date = seq(as.Date("2020-03-01"), publication_date, by = 1)), all = TRUE)

    DT[is.na(N), N := 0]

    DT[, publication_date := publication_date]

    setkey(DT, publication_date, date)

    return(DT)
}

load_fhm_icu <- function(f) {
    require(data.table)
    require(readxl)
    require(stringr)

    # Skip early reports that do not contain ICU data
    date <- as.Date(str_extract(f, "[0-9]{4}-[0-9]{2}-[0-9]{2}"))
    if (date <= as.Date("2020-04-24")) return(NULL)

    sheets <- excel_sheets(f)
    DT <- data.table((
        read_excel(path = f, sheet = grep("intensivvårdade", sheets), col_types = c("text", "numeric"))
    ))
    setnames(DT, c("date", "N"))
    DT[(tolower(date) %in% c("uppgift saknas", "uppgift saknaa", "uppgift saknas+a1")), date := NA]

    if (can_be_numeric(DT[, date])) {
        DT[, date := as.Date(as.numeric(date), origin = "1899-12-30")]
    } else {
        DT[, date := as.Date(date)]
    }

    # Ensure starting point is March 1st, and that all dates have a value
    publication_date <- get_record_date(f)
    DT <- merge(DT, data.table(date = seq(as.Date("2020-03-01"), publication_date, by = 1)), all = TRUE)
    DT[is.na(N), N := 0]

    DT[, publication_date := publication_date]

    setkey(DT, publication_date, date)

    return(DT)
}

files <- list.files(file.path("data", "FHM"), full.names = TRUE)
dts <- lapply(files, load_fhm_deaths)

deaths_dt <- rbindlist(dts)
setkey(deaths_dt, publication_date, date)

n_workdays <- Vectorize(function(a, b) {
    dates <- seq(a, b, "days")
    dates <- dates[dates != a]
    dates <- dates[!(weekdays(dates) %in% c("Saturday", "Sunday"))]
    # Drop Swedish holidays
    dates <- dates[!(dates %in% as.Date(c("2020-04-10", "2020-04-13", "2020-05-01", "2020-05-21")))]
    return(length(dates))
})

deaths_dt[!is.na(date) & publication_date > "2020-04-02", days_lag := publication_date - date]
deaths_dt[!is.na(date) & publication_date > "2020-04-02", workdays_lag := n_workdays(date, publication_date)]

deaths_dt[date == "2020-04-02" & publication_date == "2020-04-02", `:=`(workdays_lag = 0, days_lag = 0)]

# Change reports to ensure no negatives
deaths_dt[!is.na(date), n_m1 := shift(N, n = 1, type = "lag", fill = 0L), by = date]
deaths_dt[!is.na(date), n_p1 := shift(N, n = 1, type = "lead", fill = NA_integer_), by = date]
while (deaths_dt[N < n_m1, .N] > 0) {
    # cat(deaths_dt[N - n_m1 < 0, .N], ", ")
    deaths_dt[N < n_m1, N := ceiling((N + n_m1) / 2)]
    deaths_dt[N > n_p1, N := ceiling((N + n_p1) / 2)]
    deaths_dt[!is.na(date), n_m1 := shift(N, n = 1, type = "lag", fill = 0L), by = date]
    deaths_dt[!is.na(date), n_p1 := shift(N, n = 1, type = "lead", fill = NA_integer_), by = date]
}

deaths_dt[!is.na(date), n_diff := N - n_m1]
deaths_dt[, c("n_m1", "n_p1") := NULL]

setkey(deaths_dt, date, publication_date)
write_fst(deaths_dt, file.path("data", "processed", "deaths_dt.fst"))

###
# ICU stats
dts <- lapply(files, load_fhm_icu)
icu_dt <- rbindlist(dts)
setkey(icu_dt, publication_date, date)

icu_dt[!is.na(date) & publication_date > "2020-04-02", days_lag := publication_date - date]
icu_dt[!is.na(date) & publication_date > "2020-04-02", workdays_lag := n_workdays(date, publication_date)]

icu_dt[date == "2020-04-02" & publication_date == "2020-04-02", `:=`(workdays_lag = 0, days_lag = 0)]

# Change reports to ensure no negatives
icu_dt[!is.na(date), n_m1 := shift(N, n = 1, type = "lag", fill = 0L), by = date]
icu_dt[!is.na(date), n_p1 := shift(N, n = 1, type = "lead", fill = NA_integer_), by = date]
while (icu_dt[N < n_m1, .N] > 0) {
    # cat(icu_dt[N - n_m1 < 0, .N], ", ")
    icu_dt[N < n_m1, N := ceiling((N + n_m1) / 2)]
    icu_dt[N > n_p1, N := ceiling((N + n_p1) / 2)]
    icu_dt[!is.na(date), n_m1 := shift(N, n = 1, type = "lag", fill = 0L), by = date]
    icu_dt[!is.na(date), n_p1 := shift(N, n = 1, type = "lead", fill = NA_integer_), by = date]
}

icu_dt[!is.na(date), n_diff := N - n_m1]
icu_dt[, c("n_m1", "n_p1") := NULL]

write_fst(icu_dt, file.path("data", "processed", "icu_dt.fst"))

