#
# Example code illustrating plotting of DART results
#
# This is not generic, unlike most other code in this folder.


check_dart_env = new.env()
source('compare_with_gauge_environment.R', local=check_dart_env)

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
    gauge_ylims[[i]] = NA #c(-1, 1) * 2
}

## END INPUT PAR

pdf('tonga_2009_03_19_Mw7.7.pdf', width=10, height=5)
check_dart_env$compare_event_with_gauge_time_series(
    event_magnitude,
    event_hypocentre,
    event_start, 
    gauge_ids, 
    gauge_data,
    plot_durations= plot_durations,
    gauge_ylims = gauge_ylims)
dev.off()


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
event_hypocentre = c(187.43, -15.78) 
event_start = strptime('2009-09-29 17:48:11', format='%Y-%m-%d %H:%M:%S', tz='Etc/GMT')

gauge_ids = c(32401.4, 32412.4, 51425.4, 51426.4, 54401.4)

gauge_data = c(
'../../../../../DATA/TIDES/DART/dart_extract/samoa__2009_09_29_Mw8.1/samoa__2009_09_29_Mw8.1_32401.csv',
'../../../../../DATA/TIDES/DART/dart_extract/samoa__2009_09_29_Mw8.1/samoa__2009_09_29_Mw8.1_32412.csv',
'../../../../../DATA/TIDES/DART/dart_extract/samoa__2009_09_29_Mw8.1/samoa__2009_09_29_Mw8.1_51425.csv',
'../../../../../DATA/TIDES/DART/dart_extract/samoa__2009_09_29_Mw8.1/samoa__2009_09_29_Mw8.1_51426.csv',
'../../../../../DATA/TIDES/DART/dart_extract/samoa__2009_09_29_Mw8.1/samoa__2009_09_29_Mw8.1_54401.csv')

plot_durations = list()
gauge_ylims = list()
for(i in 1:length(gauge_ids)){
    plot_durations [[i]] = c(0, 3600 * 18)
    gauge_ylims[[i]] = NA #c(-1, 1) * 2
}

## END INPUT PAR

pdf('samoa_2009_09_29_Mw8.1.pdf', width=10, height=5)
check_dart_env$compare_event_with_gauge_time_series(
    event_magnitude,
    event_hypocentre,
    event_start, 
    gauge_ids, 
    gauge_data,
    plot_durations= plot_durations,
    gauge_ylims = gauge_ylims)
dev.off()
