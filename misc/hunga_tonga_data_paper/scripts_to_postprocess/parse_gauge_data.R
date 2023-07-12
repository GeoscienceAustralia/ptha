#
# Functions for reading diverse tide-gauge formats (and the single mean sea
# level pressure gauge format)
#


# This can read all the BOM mean sea level pressure data
read_pressure_gauge_BOM<-function(filename){
    x = read.csv(filename)
    time = strptime(x[,1], format='%Y-%m-%d %H:%M:%S', tz='Etc/UTC') # Already in UTC
    pressure = x[,2] # MSLPressure
    quality_flag = x[,3]
    output = data.frame(time=time, juliant=julian(time), pressure=pressure, quality_flag=quality_flag)

    k = which(!is.na(output$pressure))
    output = output[k,]

    return(output)
}

read_pressure_gauge_DES<-function(filename){
    x = read.csv(filename)
    time = -as.difftime(10, units='hours') + 
        strptime(x[,1], format='%Y-%m-%dT%H:%M', tz='Etc/UTC')
    pressure = x[,2]
    quality_flag = rep(NA, length(pressure)) # Add this for consistency with BOM pressure

    output = data.frame(time=time, juliant=julian(time), pressure=pressure, quality_flag=quality_flag)

    k = which(!is.na(output$pressure))
    output = output[k,]

    return(output)
}

# Read any pressure gauge file with the right reader, based on the filename
read_pressure_gauge_any<-function(data_file){
    if(grepl('BOM_mslp', data_file)){
        read_mslp = read_pressure_gauge_BOM
    }else if(grepl('DES_QGHL_', data_file)){
        read_mslp = read_pressure_gauge_DES
    }else{
        stop('unrecognized file type')
    }
    output = read_mslp(data_file) 
    return(output)
}


# Parse the BOM tide-gauge data
read_tide_gauge_BOM<-function(filename){
    x = read.csv(filename)
    time = strptime(x[,1], format='%Y-%m-%d %H:%M:%S', tz='Etc/UTC') # Already in UTC
    stage = x[,2]
    output = data.frame(time=time, juliant=julian(time), stage=stage)

    k = which(!is.na(output$stage))
    output = output[k,]

    return(output)
}

# This is for the TASSIE tide gauge data
read_tide_gauge_TAS<-function(filename){
    x = read.csv(filename, skip=18, header=TRUE)
    # Time is in Australian Eastern Daylight Time -- GMT+11
    # Here we convert to GMT by subtracting 11 hours
    time = strptime(x[,1], format='%d/%m/%Y %H:%M', tz='Etc/UTC') - as.difftime(11, units='hours')
    stage = x$WL_AHD
    output = data.frame(time=time, juliant=julian(time), stage=stage)

    k = which(!is.na(output$stage))
    output = output[k,]

    return(output)
}

# This is for the second round of TASSIE tide gauge data
read_tide_gauge_TAS_revised<-function(filename){
    x = read.csv(filename, skip=26, header=TRUE)
    # Time is in Australian Eastern Standard Time -- GMT+10 (changed from previous TAS data delivery)
    # Here we convert to GMT by subtracting 10 hours
    time = strptime(x[,1], format='%d/%m/%Y %H:%M', tz='Etc/UTC') - as.difftime(10, units='hours')
    stage = x$WL_AHD
    output = data.frame(time=time, juliant=julian(time), stage=stage)

    k = which(!is.na(output$stage))
    output = output[k,]

    return(output)
}

# For macquarie island
read_tide_gauge_macquarieisland<-function(filename){
    x = read.csv(filename, header=TRUE)
    # Here the timezone is already UTC
    time = strptime(x[,1], format='%d/%m/%Y %H:%M', tz='Etc/UTC') 
    stage = x[,2]
    output = data.frame(time=time, juliant=julian(time), stage=stage)
    return(output)
}

# This is for the MHL tide-gauge data
read_tide_gauge_MHL<-function(filename){
    x = read.csv(filename, skip=30, header=TRUE, na.strings='---')
    # Here the Timezone is UTC+10 -- convert to UTC
    time = strptime(paste0(x[,1], ' ', x[,2]), format='%d/%m/%Y %H:%M:%S', tz='Etc/UTC') - as.difftime(10, units='hours')
    stage = as.numeric(x[,3])
    output = data.frame(time=time, juliant=julian(time), stage=stage)

    k = which(!is.na(output$stage))
    output = output[k,]

    return(output)
}

# This is for the DES tide-gauge data
read_tide_gauge_DES<-function(filename){
    x = read.csv(filename, skip=36, header=TRUE)
    # Here the timezone is UTC+10 -- convert to UTC
    time = strptime(x[,1], format='%Y-%m-%dT%H:%M', tz='Etc/UTC') - as.difftime(10, units='hours')
    stage = as.numeric(x[,2])
    output = data.frame(time=time, juliant=julian(time), stage=stage)

    k = which(!is.na(output$stage))
    output = output[k,]

    return(output)
}

read_tide_gauge_IOC<-function(filename){
    x = read.csv(filename, header=TRUE)
    time = strptime(x[,1], format='%Y-%m-%d %H:%M:%S', tz='Etc/UTC') # Already in UTC format
    stage = as.numeric(x[,2])
    output = data.frame(time=time, juliant=julian(time), stage=stage)

    k = which(!is.na(output$stage))
    output = output[k,]
    return(output)
}

read_tide_gauge_PortAuthority<-function(filename){ 
    # Here the data is on the last line of the file
    tide_data = read.csv(filename, header=TRUE, comment.char='#', na.strings=c('NA', '---'))

    # Times are in UTC
    time_string = tide_data[,1]
    time = strptime(time_string, format='%Y-%m-%d %H:%M:%S', tz='Etc/UTC')
    stage = as.numeric(tide_data[,3])

    output = data.frame(time=time, juliant = julian(time), stage=stage)
    return(output)
}

# Try to read any tide-gauge file type
read_tide_gauge_any<-function(data_file){
    # Determine which file reader to use

    #if(grepl('BOMPorts_', data_file) | grepl('BOMQC_', data_file)){
    if(grepl('BOMPorts_', data_file)){
        read_tg = read_tide_gauge_BOM
    }else if(grepl('MHL', data_file)){
        read_tg = read_tide_gauge_MHL
    }else if(grepl('2022TsunamiTAS', data_file)){
        read_tg = read_tide_gauge_TAS
    }else if(grepl('2022_Tsunami_TAS', data_file)){
        read_tg = read_tide_gauge_TAS_revised
    }else if(grepl('Macquarie_Island', data_file)){
        read_tg = read_tide_gauge_macquarieisland
    }else if(grepl('DES_QGHL', data_file)){
        read_tg = read_tide_gauge_DES
    }else if(grepl('ioc_sealevelmonitoring', data_file)){
        read_tg = read_tide_gauge_IOC
    }else if(grepl('NSW_Port_Authority_DO_NOT_DISTRIBUTE', data_file)){
        read_tg = read_tide_gauge_PortAuthority
    }else{
        stop(paste0('Could not find reader for file ', data_file))
    }
    output = read_tg(data_file)
    return(output)
}

