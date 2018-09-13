check_dart_env = new.env()
source('compare_with_data_environment.R', local=check_dart_env)


## INPUT PAR
event_magnitude = 8.0
event_hypocentre = c(185.88, -20.19) 
event_start = strptime('2006-05-03 15:26:40', format='%Y-%m-%d %H:%M:%S', tz='Etc/GMT')
gauge_ids = 51407.4
# NOTE: Here we use '51407B', which is an 'alternative high sampling rate' data source. The typical
# data source only samples at 15min, which is too infrequent
gauge_data = '../../../../../DATA/TIDES/DART/dart_extract/tonga_2006_05_03_Mw8.0/tonga_2006_05_03_Mw8.0_51407B.csv'

plot_durations = list()
gauge_ylims = list()
for(i in 1:length(gauge_ids)){
    plot_durations [[i]] = c(0, 3600 * 24)
    gauge_ylims[[i]] = NA 
}

event_basename = 'tonga_2006_05_03_Mw8.0'
## END INPUT PAR

source('check_dart_include.R')

############################################################################
#
# Another event -- note the exact GCMT Mw is < 7.7. 
#
############################################################################


## INPUT PAR
event_magnitude = 7.7 
event_hypocentre = c(185.34, -23.05) 
event_start = strptime('2009-03-19 18:17:40', format='%Y-%m-%d %H:%M:%S', tz='Etc/GMT')

gauge_ids = 51426.4

gauge_data = '../../../../../DATA/TIDES/DART/dart_extract/tonga_2009_03_19_Mw7.7/tonga_2009_03_19_Mw7.7_51426.csv'

plot_durations = list()
gauge_ylims = list()
for(i in 1:length(gauge_ids)){
    plot_durations [[i]] = c(0, 3600 * 6)
    gauge_ylims[[i]] = NA 
}

event_basename = 'tonga_2009_03_19_Mw7.7'
## END INPUT PAR

source('check_dart_include.R')

############################################################################
#
# Another event
#
############################################################################


## INPUT PAR
event_magnitude = 8.1 
#event_hypocentre = c(187.9, -15.49) 
## NOTE: This event was outer-rise, not on our unit sources, so we deliberatly shift the 
## location to match the nearest unit source
event_hypocentre = c(360 - 172.1 - 0.6, -15.49) 
event_start = strptime('2009-09-29 17:48:11', format='%Y-%m-%d %H:%M:%S', tz='Etc/GMT')

gauge_ids = c(32401.4, 32412.4, 51425.4, 51426.4, 54401.4)

gauge_data = c(
'../../../../../DATA/TIDES/DART/dart_extract/samoa_2009_09_29_Mw8.1/samoa_2009_09_29_Mw8.1_32401.csv',
'../../../../../DATA/TIDES/DART/dart_extract/samoa_2009_09_29_Mw8.1/samoa_2009_09_29_Mw8.1_32412.csv',
'../../../../../DATA/TIDES/DART/dart_extract/samoa_2009_09_29_Mw8.1/samoa_2009_09_29_Mw8.1_51425.csv',
'../../../../../DATA/TIDES/DART/dart_extract/samoa_2009_09_29_Mw8.1/samoa_2009_09_29_Mw8.1_51426.csv',
'../../../../../DATA/TIDES/DART/dart_extract/samoa_2009_09_29_Mw8.1/samoa_2009_09_29_Mw8.1_54401.csv')

plot_durations = list()
gauge_ylims = list()
for(i in 1:length(gauge_ids)){
    plot_durations [[i]] = c(0, 3600 * 18)
    gauge_ylims[[i]] = NA #c(-1, 1) * 2
}

event_basename = 'samoa_2009_09_29_Mw8.1'
## END INPUT PAR
source('check_dart_include.R')

