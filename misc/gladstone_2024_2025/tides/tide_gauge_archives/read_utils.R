
#' Read csv files with tide gauge data and combine them into a single data frame
#' 
#' Try reading the data in the new format, if that fails, try the old format
read_data <- function(files){
    # Create data frame with data and value columns
    data <- data.frame(date=as.POSIXct(character()), value=numeric())
    errors <- list()
    for(file in files){
        print(file)
        data <- tryCatch({
            rbind(data, read_new_data(file))
        }, error=function(e){
            errors <- c(errors, e)
            data <- tryCatch({
                rbind(data, read_data_old(file))
            }, error=function(e){
                errors <- c(errors, e)
                data <- tryCatch({
                    rbind(data, read_new_data_short(file))
                }, error=function(e){
                    errors <- c(errors, e)
                    print(errors)
                    stop("Could not read data")
                })
            })
        })
    }
    check_data_read(data)
    return(data)   
}


#' Read data like this:
#' 
#' _id,010120170000  3.791
#' 1,010120170010  3.703
#' 2,010120170020  3.602
#' 43068,271020150200  0.099
# 43069,271020150210  0.049
# 43070,271020150220 -0.010
# 43071,271020150230 -0.009
# 43072,271020150240 -0.027
read_data_old <- function(file){
    data_in <- read.table(file, sep = "", header = FALSE)
    date_time_str <- sapply(strsplit(data_in$V1, ","), function(x) x[2])
    date_time <- as.POSIXct(date_time_str, format="%d%m%Y%H%M")
    value <- data_in$V2

    data <- data.frame(date=date_time, value=value)
    check_data_read
    return(data)
}

check_data_read <- function(data){
    stopifnot(!any(is.na(data$date)))
    stopifnot(!any(is.na(data$value)))
    stopifnot(nrow(data) > 0)
}

#' Read data with a 39 line header with a Date, Time, and Reading column
read_new_data <- function(file){
    data_in <- read.table(file, skip=39, header=TRUE, sep=",")
    date_time_str <- paste(data_in$Date, data_in$Time)
    date_time <- as.POSIXct(date_time_str, format="%d/%m/%Y %H:%M")

    value <- data_in$Reading
    data <- data.frame(date=date_time, value=value)
    check_data_read(data)
    return(data)
}

read_new_data_short <- function(file){
    data_in <- read.table(file, skip=0, header=TRUE, sep=",")
    date_time_str <- paste(data_in$Date, data_in$Time)
    date_time <- as.POSIXct(date_time_str, format="%d/%m/%Y %H:%M")

    value <- data_in$Reading
    data <- data.frame(date=date_time, value=value)
    check_data_read(data)
    return(data)
}

#' Read it as a csv file with two columns: names and values
read_metadata <- function(ashx_file) {
    metadata <- read.table(
        ashx_file,
        skip=31,
        sep=",",
        row.names=1,
        header=FALSE,
        nrows=16
    )
    return(metadata)
}
