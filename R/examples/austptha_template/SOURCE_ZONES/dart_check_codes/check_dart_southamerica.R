check_dart_env = new.env()
source('compare_with_data_environment.R', local=check_dart_env)

###############################################################################
#
# 2014 Mw 8.2 event
#
###############################################################################


## INPUT PAR 
event_magnitude = 8.2
event_hypocentre = c(289.23, -19.61) 
event_start = strptime('2014-04-01 23:46:47', format='%Y-%m-%d %H:%M:%S', tz='Etc/GMT')

gauge_ids = c(
    32401.4,
    32402.4,
    32412.4,
    32413.4,
    51407.4,
    51426.4,
    52406.4
    )

gauge_data = c( 
    '../../../../../DATA/TIDES/DART/dart_extract/southamerica_2014_04_01_Mw8.2/southamerica_2014_04_01_Mw8.2_32401.csv',
    '../../../../../DATA/TIDES/DART/dart_extract/southamerica_2014_04_01_Mw8.2/southamerica_2014_04_01_Mw8.2_32402.csv',
    '../../../../../DATA/TIDES/DART/dart_extract/southamerica_2014_04_01_Mw8.2/southamerica_2014_04_01_Mw8.2_32412.csv',
    '../../../../../DATA/TIDES/DART/dart_extract/southamerica_2014_04_01_Mw8.2/southamerica_2014_04_01_Mw8.2_32413.csv',
    '../../../../../DATA/TIDES/DART/dart_extract/southamerica_2014_04_01_Mw8.2/southamerica_2014_04_01_Mw8.2_51407.csv',
    '../../../../../DATA/TIDES/DART/dart_extract/southamerica_2014_04_01_Mw8.2/southamerica_2014_04_01_Mw8.2_51426.csv',
    '../../../../../DATA/TIDES/DART/dart_extract/southamerica_2014_04_01_Mw8.2/southamerica_2014_04_01_Mw8.2_52406.csv'
    )

plot_durations = list()
gauge_ylims = list()
for(i in 1:length(gauge_ids)){
    plot_durations [[i]] = c(0, 3600 * 18)
    gauge_ylims[[i]] = NA #c(-1, 1) * 2
}

event_basename = 'southamerica_2014_04_01_Mw8.2'

## END INPUT PAR

source('check_dart_include.R')

#stop()

##################################################################################
#
# An Mw 7.8 in April 2016
#
##################################################################################

event_magnitude = 7.8
event_hypocentre = c(360-79.93, 0.35)
event_start = strptime('2016-04-16 23:58:36', format='%Y-%m-%d %H:%M:%S', tz='Etc/GMT')
gauge_ids = c(32411.4, 32413.4)
gauge_data = c(
    '../../../../../DATA/TIDES/DART/dart_extract/southamerica_2016_04_16_Mw7.8/southamerica_2016_04_16_Mw7.8_32411.csv',
    '../../../../../DATA/TIDES/DART/dart_extract/southamerica_2016_04_16_Mw7.8/southamerica_2016_04_16_Mw7.8_32413.csv')
plot_durations = list()
gauge_ylims = list()
for(i in 1:length(gauge_ids)){
    plot_durations [[i]] = c(0, 3600 * 3)
    gauge_ylims[[i]] = NA
}
event_basename = 'southamerica_2016_04_16_Mw7.8'
source('check_dart_include.R')

#stop()

################################################################################
#
# A Mw 7.8 event in November 2007
#
################################################################################
event_magnitude = 7.8
event_hypocentre = c(360-69.89, -22.25) 
event_start = strptime('2007-11-14 15:40:50', format='%Y-%m-%d %H:%M:%S', tz='Etc/GMT')
gauge_ids = c(32401.4, 32412.4)
gauge_data = c(
    '../../../../../DATA/TIDES/DART/dart_extract/southamerica_2007_11_14_Mw7.8/southamerica_2007_11_14_Mw7.8_32401.csv',
    '../../../../../DATA/TIDES/DART/dart_extract/southamerica_2007_11_14_Mw7.8/southamerica_2007_11_14_Mw7.8_32412.csv')

plot_durations = list()
gauge_ylims = list()
for(i in 1:length(gauge_ids)){
    plot_durations [[i]] = c(0, 3600 * 5)
    gauge_ylims[[i]] = NA
}

event_basename = 'southamerica_2007_11_14_Mw7.8'
## END INPUT PAR

source('check_dart_include.R')

#stop()
######################################################################################
#
# 2010 Mw 8.8
#
######################################################################################


## INPUT PAR
event_magnitude = 8.8
event_hypocentre = c(287.29, -35.85) 
event_start = strptime('2010-02-27 06:34:15', format='%Y-%m-%d %H:%M:%S', tz='Etc/GMT')

gauge_ids = c(
21413.4,
32412.4,
43412.4,
46403.4,
46404.4,
46407.4,
46409.4,
46419.4,
51425.4,
51426.4,
52401.4,
52402.4,
52403.4,
52405.4,
54401.4,
55012.4)

gauge_data = c(
'../../../../../DATA/TIDES/DART/dart_extract/chile_2010_02_27_Mw8.8/chile_2010_02_27_Mw8.8_21413.csv',
'../../../../../DATA/TIDES/DART/dart_extract/chile_2010_02_27_Mw8.8/chile_2010_02_27_Mw8.8_32412.csv',
'../../../../../DATA/TIDES/DART/dart_extract/chile_2010_02_27_Mw8.8/chile_2010_02_27_Mw8.8_43412.csv',
'../../../../../DATA/TIDES/DART/dart_extract/chile_2010_02_27_Mw8.8/chile_2010_02_27_Mw8.8_46403.csv',
'../../../../../DATA/TIDES/DART/dart_extract/chile_2010_02_27_Mw8.8/chile_2010_02_27_Mw8.8_46404.csv',
'../../../../../DATA/TIDES/DART/dart_extract/chile_2010_02_27_Mw8.8/chile_2010_02_27_Mw8.8_46407.csv',
'../../../../../DATA/TIDES/DART/dart_extract/chile_2010_02_27_Mw8.8/chile_2010_02_27_Mw8.8_46409.csv',
'../../../../../DATA/TIDES/DART/dart_extract/chile_2010_02_27_Mw8.8/chile_2010_02_27_Mw8.8_46419.csv',
'../../../../../DATA/TIDES/DART/dart_extract/chile_2010_02_27_Mw8.8/chile_2010_02_27_Mw8.8_51425.csv',
'../../../../../DATA/TIDES/DART/dart_extract/chile_2010_02_27_Mw8.8/chile_2010_02_27_Mw8.8_51426.csv',
'../../../../../DATA/TIDES/DART/dart_extract/chile_2010_02_27_Mw8.8/chile_2010_02_27_Mw8.8_52401.csv',
'../../../../../DATA/TIDES/DART/dart_extract/chile_2010_02_27_Mw8.8/chile_2010_02_27_Mw8.8_52402.csv',
'../../../../../DATA/TIDES/DART/dart_extract/chile_2010_02_27_Mw8.8/chile_2010_02_27_Mw8.8_52403.csv',
'../../../../../DATA/TIDES/DART/dart_extract/chile_2010_02_27_Mw8.8/chile_2010_02_27_Mw8.8_52405.csv',
'../../../../../DATA/TIDES/DART/dart_extract/chile_2010_02_27_Mw8.8/chile_2010_02_27_Mw8.8_54401.csv',
'../../../../../DATA/TIDES/DART/dart_extract/chile_2010_02_27_Mw8.8/chile_2010_02_27_Mw8.8_55012.csv')

##/short/w85/tsunami/DATA/TIDES/DART/dart_extract/chile_2010_02_27_Mw8.8/chile_2010_02_27_Mw8.8_21413.csv
##/short/w85/tsunami/DATA/TIDES/DART/dart_extract/chile_2010_02_27_Mw8.8/chile_2010_02_27_Mw8.8_32412.csv
##/short/w85/tsunami/DATA/TIDES/DART/dart_extract/chile_2010_02_27_Mw8.8/chile_2010_02_27_Mw8.8_43412.csv
##/short/w85/tsunami/DATA/TIDES/DART/dart_extract/chile_2010_02_27_Mw8.8/chile_2010_02_27_Mw8.8_46404.csv
##/short/w85/tsunami/DATA/TIDES/DART/dart_extract/chile_2010_02_27_Mw8.8/chile_2010_02_27_Mw8.8_46407.csv
##/short/w85/tsunami/DATA/TIDES/DART/dart_extract/chile_2010_02_27_Mw8.8/chile_2010_02_27_Mw8.8_46409.csv
##/short/w85/tsunami/DATA/TIDES/DART/dart_extract/chile_2010_02_27_Mw8.8/chile_2010_02_27_Mw8.8_46419.csv
##/short/w85/tsunami/DATA/TIDES/DART/dart_extract/chile_2010_02_27_Mw8.8/chile_2010_02_27_Mw8.8_51425.csv
##/short/w85/tsunami/DATA/TIDES/DART/dart_extract/chile_2010_02_27_Mw8.8/chile_2010_02_27_Mw8.8_51426.csv
##/short/w85/tsunami/DATA/TIDES/DART/dart_extract/chile_2010_02_27_Mw8.8/chile_2010_02_27_Mw8.8_52401.csv
##/short/w85/tsunami/DATA/TIDES/DART/dart_extract/chile_2010_02_27_Mw8.8/chile_2010_02_27_Mw8.8_52402.csv
##/short/w85/tsunami/DATA/TIDES/DART/dart_extract/chile_2010_02_27_Mw8.8/chile_2010_02_27_Mw8.8_52403.csv
##/short/w85/tsunami/DATA/TIDES/DART/dart_extract/chile_2010_02_27_Mw8.8/chile_2010_02_27_Mw8.8_52405.csv
##/short/w85/tsunami/DATA/TIDES/DART/dart_extract/chile_2010_02_27_Mw8.8/chile_2010_02_27_Mw8.8_54401.csv
##/short/w85/tsunami/DATA/TIDES/DART/dart_extract/chile_2010_02_27_Mw8.8/chile_2010_02_27_Mw8.8_55012.csv



plot_durations = list()
gauge_ylims = list()
for(i in 1:length(gauge_ids)){
    plot_durations [[i]] = c(0, 3600 * 24)
    gauge_ylims[[i]] = NA #c(-1, 1) * 2
}

event_basename = 'chile_2010_02_27_Mw8.8'
## END INPUT PAR

source('check_dart_include.R')



###############################################################################
#
# Another event
#
###############################################################################


## INPUT PAR 
event_magnitude = 8.3
event_hypocentre = c(288.33, -31.57) 
event_start = strptime('2015-09-16 22:54:32', format='%Y-%m-%d %H:%M:%S', tz='Etc/GMT')

gauge_ids = c(
    #21413.4, ## < 2cm range -- remove -- because there are other 'very small' ones that we avoid as well.
    21414.4,
    21415.4,
    32402.4,
    32411.4,
    32412.4,
    43412.4,
    46403.4,
    46408.4,
    46409.4,
    46413.4,
    51407.4,
    51425.4,
    52401.4,
    52402.4,
    52403.4,
    52404.5,
    52405.4,
    52406.4
    )

gauge_data = c( 
    #'../../../../../DATA/TIDES/DART/dart_extract/southamerica_2015_09_16_Mw8.3/southamerica_2015_09_16_Mw8.3_21413.csv',
    '../../../../../DATA/TIDES/DART/dart_extract/southamerica_2015_09_16_Mw8.3/southamerica_2015_09_16_Mw8.3_21414.csv',
    '../../../../../DATA/TIDES/DART/dart_extract/southamerica_2015_09_16_Mw8.3/southamerica_2015_09_16_Mw8.3_21415.csv',
    '../../../../../DATA/TIDES/DART/dart_extract/southamerica_2015_09_16_Mw8.3/southamerica_2015_09_16_Mw8.3_32402.csv',
    '../../../../../DATA/TIDES/DART/dart_extract/southamerica_2015_09_16_Mw8.3/southamerica_2015_09_16_Mw8.3_32411.csv',
    '../../../../../DATA/TIDES/DART/dart_extract/southamerica_2015_09_16_Mw8.3/southamerica_2015_09_16_Mw8.3_32412.csv',
    '../../../../../DATA/TIDES/DART/dart_extract/southamerica_2015_09_16_Mw8.3/southamerica_2015_09_16_Mw8.3_43412.csv',
    '../../../../../DATA/TIDES/DART/dart_extract/southamerica_2015_09_16_Mw8.3/southamerica_2015_09_16_Mw8.3_46403.csv',
    '../../../../../DATA/TIDES/DART/dart_extract/southamerica_2015_09_16_Mw8.3/southamerica_2015_09_16_Mw8.3_46408.csv',
    '../../../../../DATA/TIDES/DART/dart_extract/southamerica_2015_09_16_Mw8.3/southamerica_2015_09_16_Mw8.3_46409.csv',
    '../../../../../DATA/TIDES/DART/dart_extract/southamerica_2015_09_16_Mw8.3/southamerica_2015_09_16_Mw8.3_46413.csv',
    '../../../../../DATA/TIDES/DART/dart_extract/southamerica_2015_09_16_Mw8.3/southamerica_2015_09_16_Mw8.3_51407.csv',
    '../../../../../DATA/TIDES/DART/dart_extract/southamerica_2015_09_16_Mw8.3/southamerica_2015_09_16_Mw8.3_51425.csv',
    '../../../../../DATA/TIDES/DART/dart_extract/southamerica_2015_09_16_Mw8.3/southamerica_2015_09_16_Mw8.3_52401.csv',
    '../../../../../DATA/TIDES/DART/dart_extract/southamerica_2015_09_16_Mw8.3/southamerica_2015_09_16_Mw8.3_52402.csv',
    '../../../../../DATA/TIDES/DART/dart_extract/southamerica_2015_09_16_Mw8.3/southamerica_2015_09_16_Mw8.3_52403.csv',
    '../../../../../DATA/TIDES/DART/dart_extract/southamerica_2015_09_16_Mw8.3/southamerica_2015_09_16_Mw8.3_52404.csv',
    '../../../../../DATA/TIDES/DART/dart_extract/southamerica_2015_09_16_Mw8.3/southamerica_2015_09_16_Mw8.3_52405.csv',
    '../../../../../DATA/TIDES/DART/dart_extract/southamerica_2015_09_16_Mw8.3/southamerica_2015_09_16_Mw8.3_52406.csv'
    )

plot_durations = list()
gauge_ylims = list()
for(i in 1:length(gauge_ids)){
    plot_durations [[i]] = c(0, 3600 * 30)
    gauge_ylims[[i]] = NA #c(-1, 1) * 2
}

event_basename = 'southamerica_2015_09_16_Mw8.3'

## END INPUT PAR

source('check_dart_include.R')



###############################################################################
#
# Another event
#
###############################################################################


## INPUT PAR 
event_magnitude = 8.0
event_hypocentre = c(283.4, -13.40) 
event_start = strptime('2007-08-15 23:40:57', format='%Y-%m-%d %H:%M:%S', tz='Etc/GMT')

gauge_ids = c(32401.4, 32411.4, 43412.4)
gauge_data = c(
    '../../../../../DATA/TIDES/DART/dart_extract/southamerica_2007_08_15_Mw8.0/southamerica_2007_08_15_Mw8.0_32401.csv',
    '../../../../../DATA/TIDES/DART/dart_extract/southamerica_2007_08_15_Mw8.0/southamerica_2007_08_15_Mw8.0_32411.csv',
    '../../../../../DATA/TIDES/DART/dart_extract/southamerica_2007_08_15_Mw8.0/southamerica_2007_08_15_Mw8.0_43412.csv')
plot_durations = list()
gauge_ylims = list()
for(i in 1:length(gauge_ids)){
    plot_durations [[i]] = c(0, 3600 * 30)
    gauge_ylims[[i]] = NA #c(-1, 1) * 2
}
event_basename = 'southamerica_2007_08_15_Mw80'

## END INPUT PAR

source('check_dart_include.R')



