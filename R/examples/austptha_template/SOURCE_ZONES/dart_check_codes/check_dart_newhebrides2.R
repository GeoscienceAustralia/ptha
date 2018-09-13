check_dart_env = new.env()
source('compare_with_data_environment.R', local=check_dart_env)

## INPUT PAR
event_magnitude = 7.8 
event_hypocentre = c(166.38, -12.52) 
event_start = strptime('2009-10-07 22:18:51', format='%Y-%m-%d %H:%M:%S', tz='Etc/GMT')

gauge_ids = 51425.4

gauge_data = '../../../../../DATA/TIDES/DART/dart_extract/vanuatu_north_2009_10_07_Mw7.8/vanuatu_north_2009_10_07_Mw7.8_51425.csv'

plot_durations = list()
gauge_ylims = list()
for(i in 1:length(gauge_ids)){
    plot_durations [[i]] = c(0, 3600 * 6)
    gauge_ylims[[i]] = NA #c(-1, 1) * 2
}

event_basename = 'vanuatu_north_2009_10_07_Mw7.8'
## END INPUT PAR
source('check_dart_include.R')

#
#
# Another event
#
#

# INPUT PAR
event_magnitude = 7.9 
event_hypocentre = c(165.11, -10.8) 
event_start = strptime('2013-02-06 01:12:25', format='%Y-%m-%d %H:%M:%S', tz='Etc/GMT')

gauge_ids = c(
    51425.4,
    52401.4,
    52402.4,
    52406.4,
    55012.4)

gauge_data = c(
    '../../../../../DATA/TIDES/DART/dart_extract/northNewHebrides_2013_02_06_Mw7.9/northNewHebrides_2013_02_06_Mw7.9_51425.csv',
    '../../../../../DATA/TIDES/DART/dart_extract/northNewHebrides_2013_02_06_Mw7.9/northNewHebrides_2013_02_06_Mw7.9_52401.csv',
    '../../../../../DATA/TIDES/DART/dart_extract/northNewHebrides_2013_02_06_Mw7.9/northNewHebrides_2013_02_06_Mw7.9_52402.csv',
    '../../../../../DATA/TIDES/DART/dart_extract/northNewHebrides_2013_02_06_Mw7.9/northNewHebrides_2013_02_06_Mw7.9_52406.csv',
    '../../../../../DATA/TIDES/DART/dart_extract/northNewHebrides_2013_02_06_Mw7.9/northNewHebrides_2013_02_06_Mw7.9_55012.csv')


plot_durations = list()
gauge_ylims = list()
for(i in 1:length(gauge_ids)){
    plot_durations [[i]] = c(0, 3600 * 6)
    gauge_ylims[[i]] = NA #c(-1, 1) * 2
}

event_basename = 'northNewHebrides_2013_02_06_Mw7.9'
## END INPUT PAR
source('check_dart_include.R')

