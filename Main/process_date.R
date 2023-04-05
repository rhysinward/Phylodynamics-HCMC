process_date <- function(x){
date <- ifelse(nchar(x$Collection_Date) == 4, paste(x$Collection_Date, "06-15", sep = "-"),
                                                    ifelse(nchar(x$Collection_Date) == 7, paste(x$Collection_Date, 15, sep = "-"), x$Collection_Date))
x$Date <- as.Date(parse_date_time(date, orders = c('mdy','dmy','myd','y','my','m','ymd','ym')))

x <- dplyr :: select(x,-c('Collection_Date'))
}
