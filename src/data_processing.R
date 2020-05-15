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
    require(readxl)
    sheets <- excel_sheets(f)
    return(as.Date(sub("^FOHM ", "", sheets[length(sheets)]), format="%d %b %Y"))
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

    DT[date == "Uppgift saknas" | date == "uppgift saknas", date := NA]

    if (can_be_numeric(DT[, date])) {
        DT[, date := as.Date(as.numeric(date), origin = "1899-12-30")]
    } else {
        DT[, date := as.Date(date)]
    }

    DT[is.na(N), N := 0]

    DT[, publication_date := get_record_date(f)]

    return(DT)
}

files <- list.files(file.path("data", "FHM"), full.names = TRUE)
dts <- lapply(files, load_fhm_deaths)

death_dt <- rbindlist(dts)
setkey(death_dt, publication_date, date)

death_dt[!is.na(date) & publication_date > "2020-04-02", days_since_publication := publication_date - date]
death_dt[date == "2020-04-02" & publication_date == "2020-04-02", days_since_publication := 0]

death_dt[!is.na(date), paste0("n_m", 1) := shift(N, n = 1, type = "lag", fill = 0L), by = date]
death_dt[!is.na(date), n_diff := N - n_m1]
death_dt[!is.na(date) & n_m1 > 0 & !is.na(n_m1), n_diff_pct := N/n_m1 - 1]
death_dt[!is.na(date) & n_m1 == 0 & N == 0, n_diff_pct := 0]
death_dt[, n_m1 := NULL]

# If no death reported on publication date
for (i in seq_along(unique(death_dt$publication_date))) {
    pub <- unique(death_dt$publication_date)[i]
    if (death_dt[date == publication_date & publication_date == pub, .N] == 0) {
        death_dt <- rbind(death_dt,
                          data.table(date = pub, N = 0,
                                     publication_date = pub,
                                     days_since_publication = as.difftime(0, units = "days"),
                                     n_diff = 0, n_diff_pct = 0))
    }
}

write_fst(death_dt, file.path("data", "processed", "deaths_dt.fst"))
