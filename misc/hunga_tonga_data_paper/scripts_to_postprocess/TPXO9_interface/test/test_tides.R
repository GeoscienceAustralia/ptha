######################################################################################################
#' Read data from a csv file containing the gauge info
#read_csv_gauge<-function(gauge_file){
#
#    JULIAN_TIME_ORIGIN = strptime("1970-01-01 00:00:00", format='%Y-%m-%d %H:%M:%S', tz = "Etc/GMT-10")
#    tomaree = read.csv(gauge_file, skip=17,header=T,
#        na.strings=c('NA', 'n/a'), stringsAsFactors=FALSE)
#    tomaree_time = strptime(paste(tomaree[,1], tomaree[,2]), 
#        format='%d/%m/%Y %H:%M:%S', tz='Etc/GMT-10')
#    tomaree_julian_time = julian(tomaree_time, 
#        origin = JULIAN_TIME_ORIGIN)
#
#    output = data.frame(time=tomaree_time, 
#                        julian_time = as.numeric(tomaree_julian_time),
#                        tide=tomaree[,3], status=tomaree[,4])
#    return(output)
#}
#
#coffsh = read_csv_gauge('../../../../../DATA/NSW_Tidal_Gauges/Nov2014_PubWeb_Data_to_NTC_Harbour/Coffs HarbourPW.csv')
#t1 = strptime('2001-10-01 00:00:00', format='%Y-%m-%d %H:%M:%S', tz='Etc/GMT-10')
#t2 = strptime('2001-11-01 00:00:00', format='%Y-%m-%d %H:%M:%S', tz='Etc/GMT-10')
#ci = which(coffsh$time >= t1 & coffsh$time <= t2)
#saveRDS(coffsh[ci,], 'coffs_data_test.RDS')
coffsh = readRDS('coffs_data_test.RDS')

########################################################################################################
# Get some predictions
source('../predict_tide.R', chdir=TRUE)

site_name = 'Coffs_harbour' # (No spaces in name)
# Site coordinates as decimal degrees (longitude, latitude)
site_coordinates = c(153 + 48/60 + 45.82/(60*60), -(30 + 18/60 + 10.29/(60*60)))

# Timezone MUST be GMT / UTC. Note that -10 means 'ahead 10 hours', 
start_time = coffsh$time[1]
end_time   = coffsh$time[length(coffsh$time)]
time_interval = '15 min'

coffs_pred = get_tidal_prediction(site_name, site_coordinates, 
    start_time, end_time, time_interval)

## Compare

# Times should be identical
stopifnot(all(coffsh$time == coffs_pred[,1]))

# Tidal residual
tidal_residual = coffsh$tide - coffs_pred[,2]
mean_offset = mean(coffsh$tide, na.rm=T)

pdf('tidal_test_plot.pdf', width=10, height=7)
plot(coffsh$time, coffsh$tide - mean_offset, t='l')
points(coffs_pred$time, coffs_pred$tide, t='l', col='red')
points(coffs_pred$time, coffs_pred$tide - (coffsh$tide - mean_offset), t='l', col='green')
legend('topright', c('Data - MSL', 'Model', 'Error'), col=c('black', 'red', 'green'), 
    lty=c(1,1,1))
dev.off()

# Test that the residual standard deviation is < 7 cm
if(sd(coffs_pred$tide + mean_offset - coffsh$tide) < 0.07){
    print('PASS')
}else{
    print('FAIL')
}

