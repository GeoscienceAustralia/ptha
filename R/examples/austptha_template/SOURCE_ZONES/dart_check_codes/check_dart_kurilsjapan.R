check_dart_env = new.env()
source('compare_with_data_environment.R', local=check_dart_env)

#
# 2006 Kuril Islands Event
#
## INPUT PAR
event_magnitude = 8.3
event_hypocentre = c(153.29, 46.57) 
event_start = strptime('2006-11-15 11:14:17', format='%Y-%m-%d %H:%M:%S', tz='Etc/GMT')

gauge_ids = c(
    21414.4,
    32401.4,
    46402.4,
    46403.4,
    46404.4,
    46408.4,
    46409.4,
    46410.4,
    46411.4,
    46412.4,
    46413.4)

gauge_data = paste0('../../../../../DATA/TIDES/DART/dart_extract/kuril_2006_11_15_Mw8.3/kuril_2006_11_15_Mw8.3_', 
    floor(gauge_ids), '.csv')

plot_durations = list()
gauge_ylims = list()
for(i in 1:length(gauge_ids)){
    plot_durations [[i]] = c(0, 3600 * 24)
    gauge_ylims[[i]] = NA #c(-1, 1) * 2
}

event_basename = 'kuril_2006_11_15_Mw8.3'

## END INPUT PAR

source('check_dart_include.R')

#stop()
#
# 2011 Tohoku Tsunami!
#

## INPUT PAR
event_magnitude = 9.1
event_hypocentre = c(142.37, 38.32) 
event_start = strptime('2011-03-11 05:46:23', format='%Y-%m-%d %H:%M:%S', tz='Etc/GMT')

gauge_ids = c(
21401.4,
21413.4,
21414.4,
21415.4,
21418.4,
21419.4,
32401.4,
32411.4,
32412.4,
32413.4,
43412.4,
43413.4,
46402.4,
46403.4,
46404.4,
46408.4,
46409.4,
46410.4,
46411.4,
46412.4,
51407.4,
51425.4,
52402.4,
52403.4,
52405.4,
52406.4,
55012.4,
55023.4)

gauge_data = c(
'../../../../../DATA/TIDES/DART/dart_extract/tohoku_2011_03_11_Mw9.1/tohoku_2011_03_11_Mw9.1_21401.csv',
'../../../../../DATA/TIDES/DART/dart_extract/tohoku_2011_03_11_Mw9.1/tohoku_2011_03_11_Mw9.1_21413.csv',
'../../../../../DATA/TIDES/DART/dart_extract/tohoku_2011_03_11_Mw9.1/tohoku_2011_03_11_Mw9.1_21414.csv',
'../../../../../DATA/TIDES/DART/dart_extract/tohoku_2011_03_11_Mw9.1/tohoku_2011_03_11_Mw9.1_21415.csv',
'../../../../../DATA/TIDES/DART/dart_extract/tohoku_2011_03_11_Mw9.1/tohoku_2011_03_11_Mw9.1_21418.csv',
'../../../../../DATA/TIDES/DART/dart_extract/tohoku_2011_03_11_Mw9.1/tohoku_2011_03_11_Mw9.1_21419.csv',
'../../../../../DATA/TIDES/DART/dart_extract/tohoku_2011_03_11_Mw9.1/tohoku_2011_03_11_Mw9.1_32401.csv',
'../../../../../DATA/TIDES/DART/dart_extract/tohoku_2011_03_11_Mw9.1/tohoku_2011_03_11_Mw9.1_32411.csv',
'../../../../../DATA/TIDES/DART/dart_extract/tohoku_2011_03_11_Mw9.1/tohoku_2011_03_11_Mw9.1_32412.csv',
'../../../../../DATA/TIDES/DART/dart_extract/tohoku_2011_03_11_Mw9.1/tohoku_2011_03_11_Mw9.1_32413.csv',
'../../../../../DATA/TIDES/DART/dart_extract/tohoku_2011_03_11_Mw9.1/tohoku_2011_03_11_Mw9.1_43412.csv',
'../../../../../DATA/TIDES/DART/dart_extract/tohoku_2011_03_11_Mw9.1/tohoku_2011_03_11_Mw9.1_43413.csv',
'../../../../../DATA/TIDES/DART/dart_extract/tohoku_2011_03_11_Mw9.1/tohoku_2011_03_11_Mw9.1_46402.csv',
'../../../../../DATA/TIDES/DART/dart_extract/tohoku_2011_03_11_Mw9.1/tohoku_2011_03_11_Mw9.1_46403.csv',
'../../../../../DATA/TIDES/DART/dart_extract/tohoku_2011_03_11_Mw9.1/tohoku_2011_03_11_Mw9.1_46404.csv',
'../../../../../DATA/TIDES/DART/dart_extract/tohoku_2011_03_11_Mw9.1/tohoku_2011_03_11_Mw9.1_46408.csv',
'../../../../../DATA/TIDES/DART/dart_extract/tohoku_2011_03_11_Mw9.1/tohoku_2011_03_11_Mw9.1_46409.csv',
'../../../../../DATA/TIDES/DART/dart_extract/tohoku_2011_03_11_Mw9.1/tohoku_2011_03_11_Mw9.1_46410.csv',
'../../../../../DATA/TIDES/DART/dart_extract/tohoku_2011_03_11_Mw9.1/tohoku_2011_03_11_Mw9.1_46411.csv',
'../../../../../DATA/TIDES/DART/dart_extract/tohoku_2011_03_11_Mw9.1/tohoku_2011_03_11_Mw9.1_46412.csv',
'../../../../../DATA/TIDES/DART/dart_extract/tohoku_2011_03_11_Mw9.1/tohoku_2011_03_11_Mw9.1_51407.csv',
'../../../../../DATA/TIDES/DART/dart_extract/tohoku_2011_03_11_Mw9.1/tohoku_2011_03_11_Mw9.1_51425.csv',
'../../../../../DATA/TIDES/DART/dart_extract/tohoku_2011_03_11_Mw9.1/tohoku_2011_03_11_Mw9.1_52402.csv',
'../../../../../DATA/TIDES/DART/dart_extract/tohoku_2011_03_11_Mw9.1/tohoku_2011_03_11_Mw9.1_52403.csv',
'../../../../../DATA/TIDES/DART/dart_extract/tohoku_2011_03_11_Mw9.1/tohoku_2011_03_11_Mw9.1_52405.csv',
'../../../../../DATA/TIDES/DART/dart_extract/tohoku_2011_03_11_Mw9.1/tohoku_2011_03_11_Mw9.1_52406.csv',
'../../../../../DATA/TIDES/DART/dart_extract/tohoku_2011_03_11_Mw9.1/tohoku_2011_03_11_Mw9.1_55012.csv',
'../../../../../DATA/TIDES/DART/dart_extract/tohoku_2011_03_11_Mw9.1/tohoku_2011_03_11_Mw9.1_55023.csv')

plot_durations = list()
gauge_ylims = list()
for(i in 1:length(gauge_ids)){
    plot_durations [[i]] = c(0, 3600 * 24)
    gauge_ylims[[i]] = NA #c(-1, 1) * 2
}

event_basename = 'tohoku_2011_03_11_Mw9.1'

## END INPUT PAR

source('check_dart_include.R')
