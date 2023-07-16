## Crude R interface to use of predict_tide
## Uses system calls, but takes care of creating input files, calling, etc.

## The actual case specific code is not here -- but these routines can be used to do the work

# Get key variables for the installed TPXO 'predict_tide' program
#source('OTPS_directory_name.R', local=TRUE)
source('OTPS_directory_name_TPXO9.R', local=TRUE)

###################################################################################
#
# BELOW HERE THINGS SHOULDN'T BE CHANGED
#
###################################################################################

#' Make the lat_lon_time file for predict_tide
.make_lat_lon_time_file<-function(site_coordinates, prediction_times, 
    output_file = NA, verbose=TRUE){

    if(verbose) print('Making lat_lon_time_file...')
   
    lp = length(prediction_times)
    time_info = gsub("-", " ", prediction_times)
    time_info = gsub(":", " ", time_info)
    output_data = paste(
        rep(site_coordinates[2], len=lp),
        rep(site_coordinates[1], len=lp),
        time_info,
        sep=" ")

    if(!is.na(output_file)){
        cat(output_data, file=output_file, sep="\n")
    }

    return()
}


#' Make the setup.inp file for predict_tide
.make_setup.inp<-function(lat_lon_time_file, output_tide_file, setup_file, 
    verbose=TRUE){
    if(verbose) print('Making setup.inp ...')
    # First line = model data source
    file_lines = .TIDAL_MODEL_CONTROL_FILE
    # Second = file with lat/lon/time
    file_lines = c(file_lines, basename(lat_lon_time_file))
    # Third = What to get
    file_lines = c(file_lines, 'z')
    # 4th/5th = blank
    file_lines = c(file_lines, '', '')
    # 5th = height datum
    file_lines = c(file_lines, 'oce')
    # 6th = option which we must set to zero
    file_lines = c(file_lines, '0')
    # 7th = output file
    file_lines = c(file_lines, basename(output_tide_file))

    # Write the output
    cat(file_lines, file=setup_file, sep="\n")
}

#' Make the command line call
.run_predict_tide<-function(setup_file, OTPS_directory, verbose=TRUE){

    if(verbose) print('Calling predict_tide ... (its output printed below) ...')
    # Go back to the current directory when this function finishes
    mydir = getwd()
    on.exit(setwd(mydir))

    setwd(OTPS_directory)

    # Run the predict_tide program
    system_command = paste0('./predict_tide <', basename(setup_file))
    system(system_command)

}

#' Get rid of temp files used to run predict_tide
.cleanup_intermediate_files<-function(lat_lon_time_file, output_tide_file, 
    setup_file){
    file.remove(lat_lon_time_file)
    file.remove(setup_file)
    file.remove(output_tide_file)
}


#' Determine if a year is a leap year
.is_leap_year<-function(year){

    (year%%4 == 0)*( (year%%100 != 0) + (year%%100 == 0)*(year%%400 == 0))
}

.TEN_THOUSAND_YEARS_DAYS = sum(.is_leap_year(0:9999))*366 + sum(!.is_leap_year(0:9999))*365

#' Timezone conversion + deal with years > 9999
#'
#' @param start_time time to convert (strptime object)
#' @param tz Timezone (e.g. 'Etc/GMT-10')
#' @return Strptime object with start_time converted to the new timezone
#'
.convert_timezone<-function(start_time, tz){
    
        time_format_string = '%Y-%m-%d %H:%M:%S'

        year = as.numeric(format(start_time, '%Y'))

        if(year >= 1e+04){

            extra_ten_thousands = floor(year/1e+04)

            extra_days = as.difftime(extra_ten_thousands*.TEN_THOUSAND_YEARS_DAYS, units='days')

            start_time = start_time - extra_days

        }else{

            extra_days = as.difftime(0, units='days')
        }

        start_time_char = format(as.POSIXct(start_time), tz=tz, 
            format=time_format_string)

        new_start_time = strptime(start_time_char, format=time_format_string,
            tz=tz)

        new_start_time = new_start_time + extra_days

        return(new_start_time)
}

###############################################################################
#
# Main program below here
#
###############################################################################

#' Get tides at a particular site using TPXO
#'
#' @param site_name Name for the site (just used for files)
#' @param site_coordinates vector of length 2 giving lon,lat in decimal degrees
#' @param start_time start time (strptime object). Must have correct timezone
#' @param end_time end time (strptime object). Must have correct timezone
#' @param time_interval. Time interval accepted by 'seq' (e.g. '1 hour' or '15 min' or '5 days' or '30 sec')
#' @param OPTS_directory location of the TPXO OPTS directory
#' @return data.frame with time and tidal level at the site coordinates
#' 
get_tidal_prediction<-function(site_name, site_coordinates, 
    start_time, end_time, time_interval, OTPS_directory=.OTPS_directory){

    stopifnot(length(start_time) == length(end_time))

    # Key filenames
    lat_lon_time_file = paste0(OTPS_directory, '/' , 'lat_lon_time_R_', 
        site_name)
    output_tide_file = paste0(OTPS_directory, '/' , 'predictions_', site_name,
        '.out')
    setup_file = paste0(OTPS_directory, '/', 'setup_', site_name, '.inp')


    prediction_times = vector(mode='list', len = length(start_time))
    gmt_prediction_times = vector(mode='list', len = length(start_time))
    tz = attr(start_time[1], 'tzone')[1] 
    format_string = '%Y-%m-%d %H:%M:%S'
    GMT_tz = 'Etc/GMT'
    apply_timezone_conversion = (tz != GMT_tz)
    for(i in 1:length(start_time)){

        if(i%%1e+04 == 0) print(i)

        if(apply_timezone_conversion){
            # Compute times to get output
            pred_seq = seq(start_time[i], end_time[i], by=time_interval)

            # Store the 'time-difference' between the desired prediction time,
            # and the FIRST start_time (which is simply used as a datum)
            # This is done to work-around some problems with converting to a vector
            # (without care, R converts the times to numeric and we lose key information)
            prediction_times[[i]] = (pred_seq - start_time[1])

            lps = length(pred_seq)

            # Change timezone
            gmt_start_time = .convert_timezone(pred_seq[1], tz=GMT_tz)
            gmt_end_time = .convert_timezone(pred_seq[lps], tz=GMT_tz)

            gmt_prediction_times[[i]] = strftime(seq(gmt_start_time, gmt_end_time, len=lps), 
                format=format_string, tz=GMT_tz)

            stopifnot(length(prediction_times[[i]]) == length(gmt_prediction_times[[i]]))
        }else{
            # Compute times to get output
            pred_seq = seq(start_time[i], end_time[i], by=time_interval)
            # Store the 'time-difference' between the desired prediction time,
            # and the FIRST start_time (which is simply used as a datum)
            # This is done to work-around some problems with converting to a vector
            # (without care, R converts the times to numeric and we lose key information)
            prediction_times[[i]] = pred_seq - start_time[1]
        }
    }

    print('Unpack 1 ...')
    # Pack the prediction times into a single vector. This is surprisingly
    # hard to do without coercian to numeric, which causes errors
    pt_lengths = unlist(lapply(prediction_times, length))
    prediction_times2 = rep(start_time[1], length = sum(pt_lengths))
    start_ind = 1
    for(i in 1:length(start_time)){
        end_ind = pt_lengths[i]
        prediction_times2[start_ind:end_ind] = prediction_times[[i]] + start_time[1]
        start_ind = end_ind + 1
    }
    prediction_times = prediction_times2

    print('Unpack 2 ...')
    if(tz == GMT_tz){
        gmt_prediction_times = strftime(prediction_times, format=format_string, tz=GMT_tz)
    }else{
        gmt_prediction_times = unlist(gmt_prediction_times) 
    }

    # Make the lat_lon_time file
    .make_lat_lon_time_file(site_coordinates, prediction_times=gmt_prediction_times, 
        output_file=lat_lon_time_file)

    # Make setup.inp
    .make_setup.inp(lat_lon_time_file, output_tide_file, setup_file)

    # Call predict_tide
    .run_predict_tide(setup_file, OTPS_directory)

    # Read the output of predict_tide
    output_text = read.table(output_tide_file, sep="", colClasses='character', 
        stringsAsFactors=FALSE, skip=6)

    output_data = data.frame(time=prediction_times, tide=as.numeric(output_text[,5]))

    .cleanup_intermediate_files(lat_lon_time_file, setup_file, output_tide_file)

    return(output_data)
}


