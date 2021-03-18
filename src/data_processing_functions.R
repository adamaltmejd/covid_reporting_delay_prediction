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

    DT[, publication_date := get_record_date(f)]

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
        read_excel(path = f, sheet = grep("intensivvårdade", sheets), col_types = c("text", "numeric"), progress = FALSE)
    ))
    setnames(DT, c("date", "N"))
    DT[(tolower(date) %in% c("uppgift saknas", "uppgift saknaa", "uppgift saknas+a1")), date := NA]

    if (can_be_numeric(DT[, date])) {
        DT[, date := as.Date(as.numeric(date), origin = "1899-12-30")]
    } else {
        DT[, date := as.Date(date)]
    }

    DT[, publication_date := get_record_date(f)]

    setkey(DT, publication_date, date)

    return(DT)
}

prepare_dt <- function(DT) {
    # Remove remports that have not been assigned a date
    DT <- DT[!is.na(date) & !is.na(publication_date)]

    DT[, publication_date := as.Date(publication_date)]
    DT[, date := as.Date(date)]

    # Ensure complete time series
    all_dates <- CJ(date = seq.Date(DT[, min(date)], DT[, max(publication_date)], by = 1),
                    publication_date = seq.Date(DT[, min(publication_date)], DT[, max(publication_date)], by = 1))
    all_dates <- all_dates[date <= publication_date]
    DT <- merge(DT, all_dates, by = names(all_dates), all = TRUE)

    # Fill dates with no publication with the same value as the day before
    DT[, N := nafill(N, type = "locf"), by = .(date)]

    # Set missing to zero (happens when a date starts on a day without a report)
    DT[is.na(N), N := 0]

    # Calculate lags
    n_workdays <- Vectorize(function(a, b) {
        dates <- seq(a, b, "days")
        dates <- dates[dates != a]
        dates <- dates[!(weekdays(dates) %in% c("Saturday", "Sunday"))]
        # Drop Swedish holidays
        dates <- dates[!(dates %in% as.Date(c("2020-04-10", "2020-04-13", "2020-05-01", "2020-05-21")))]
        return(length(dates))
    })

    DT[, days_lag := publication_date - date]
    DT[, workdays_lag := n_workdays(date, publication_date)]

    # Change reports to ensure no negatives
    DT[!is.na(date), n_m1 := shift(N, n = 1, type = "lag", fill = 0L), by = date]
    DT[!is.na(date), n_p1 := shift(N, n = 1, type = "lead", fill = NA_integer_), by = date]
    while (DT[N < n_m1, .N] > 0) {
        # cat(DT[N - n_m1 < 0, .N], ", ")
        DT[N < n_m1, N := ceiling((N + n_m1) / 2)]
        DT[N > n_p1, N := ceiling((N + n_p1) / 2)]
        DT[!is.na(date), n_m1 := shift(N, n = 1, type = "lag", fill = 0L), by = date]
        DT[!is.na(date), n_p1 := shift(N, n = 1, type = "lead", fill = NA_integer_), by = date]
    }

    DT[!is.na(date), n_diff := N - n_m1]
    DT[, c("n_m1", "n_p1") := NULL]

    # Set n_diff to zero when there was no report
    DT[is.na(n_diff), n_diff := 0]

    setkey(DT, date, publication_date)

    return(DT)
}
