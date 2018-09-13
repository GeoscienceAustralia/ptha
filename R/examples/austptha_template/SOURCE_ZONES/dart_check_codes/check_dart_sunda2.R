check_dart_env = new.env()
source('compare_with_data_environment.R', local=check_dart_env)

## INPUT PAR
event_magnitude = 8.5
event_hypocentre = c(101.37, -4.44) 
event_start = strptime('2007-09-12 11:10:26', format='%Y-%m-%d %H:%M:%S', tz='Etc/GMT')

gauge_ids = 23401.4

gauge_data = '../../../../../DATA/TIDES/DART/dart_extract/sumatra_2007_09_12_Mw8.5/sumatra_2007_09_12_Mw8.5_23401.csv'

plot_durations = list()
gauge_ylims = list()
for(i in 1:length(gauge_ids)){
    plot_durations [[i]] = c(0, 3600 * 24)
    gauge_ylims[[i]] = NA #c(-1, 1) * 2
}

event_basename = 'sumatra_2007_09_12_Mw8.5'
## END INPUT PAR
source('check_dart_include.R')

##############################################################################
#
# Another Event
#
##############################################################################

## INPUT PAR
event_magnitude = 7.8
event_hypocentre = c(97.05, 2.38) 
event_start = strptime('2010-04-06 22:15:01', format='%Y-%m-%d %H:%M:%S', tz='Etc/GMT')

#gauge_ids = c(23401.4, 56001.4)  # On reflection, 56001 is at 'noise level'
gauge_ids = c(23401.4) #, 56001.4) 

gauge_data = c(
    '../../../../../DATA/TIDES/DART/dart_extract/sumatra_2010_04_06_Mw7.8/sumatra_2010_04_06_Mw7.8_23401.csv')#,
    #'../../../../../DATA/TIDES/DART/dart_extract/sumatra_2010_04_06_Mw7.8/sumatra_2010_04_06_Mw7.8_56001.csv')

plot_durations = list()
gauge_ylims = list()
for(i in 1:length(gauge_ids)){
    plot_durations [[i]] = c(0, 3600 * 24)
    gauge_ylims[[i]] = NA #c(-1, 1) * 2
}

event_basename = 'sumatra_2010_04_06_Mw7.8'
## END INPUT PAR
source('check_dart_include.R')
 
##############################################################################
#
# Another Event
#
##############################################################################

## INPUT PAR
event_magnitude = 7.9
event_hypocentre = c(100.08, -3.49) 
event_start = strptime('2010-10-25 14:42:22', format='%Y-%m-%d %H:%M:%S', tz='Etc/GMT')

gauge_ids = 56001.4

gauge_data = '../../../../../DATA/TIDES/DART/dart_extract/mentawai_2010_10_25_Mw7.9/mentawai_2010_10_25_Mw7.9_56001.csv'

plot_durations = list()
gauge_ylims = list()
for(i in 1:length(gauge_ids)){
    plot_durations [[i]] = c(0, 3600 * 24)
    gauge_ylims[[i]] = NA #c(-1, 1) * 2
}
event_basename = 'mentawai_2010_10_25_Mw7.9'
## END INPUT PAR
source('check_dart_include.R')

