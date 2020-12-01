# This code provides a consistent interface to the different tide-gauge records. It is used in out 
# plotting scripts [e.g. ../../swals/plots] where the code is greatly
# simplified by having a consistent interface.

# Setup the data filepaths in a way that the function can be called from any directory on the machine
gauge_files = list(
    'Tonga2006'  = list(
        file=normalizePath('Tonga2006_BOM/nukualofa_2006_tsunami_detided.csv'),
        start_time=strptime('2006/05/03 15:26:40', format='%Y/%m/%d %H:%M:%S', tz='Etc/GMT')),
    'Tohoku2011' = list(
        file=normalizePath('Tohoku2011_IOC_Sealevel/nukualofa_tohoku_tsunami_detided.csv'),
        start_time=strptime('2011/03/11 05:46:23', format='%Y/%m/%d %H:%M:%S', tz='Etc/GMT')),
    'Chile2010' = list(
        file=normalizePath('BOM_2009_2010/nukualofa_chile2010_detided.csv'),
        start_time=strptime('2010/02/27 06:34:15', format='%Y/%m/%d %H:%M:%S', tz='Etc/GMT')),
    'Chile2015' = list(
        file=normalizePath('BOM_2014_2020/nukualofa_chile2015_detided.csv'),
        start_time=strptime('2015/09/16 22:54:32', format='%Y/%m/%d %H:%M:%S', tz='Etc/GMT'))
    )

# Get the data in a consistent format, given an event_name [corresponding to a name in the gauge_files list]
get_gauge_data<-function(event_name){

    target_file = gauge_files[[event_name]]$file
    if(!file.exists(target_file)){
        stop(paste0('Could not find file matching event_name: ', event_name))
    }

    data = read.csv(target_file)
    data$time = strptime(data$time, format='%Y-%m-%d %H:%M:%S', tz='Etc/GMT')

    # Include the coordinate for output
    event_data = list()
    event_data$obs = data
    # Coordinate of Nuku'alofa tide-gage, based on BOM metadata
    event_data$coord = c(360 -175.1815, -21.1380)
    event_data$gauge_name = "nukualofa"
    event_data$start_time = gauge_files[[event_name]]$start_time
    return(event_data)
}
