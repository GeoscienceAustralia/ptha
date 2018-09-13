check_dart_env = new.env()
source('compare_with_data_environment.R', local=check_dart_env)

## INPUT PAR
event_magnitude = 8.1
event_hypocentre = c(157.04, -8.46)
event_start = strptime('2007-04-01 20:39:56', format='%Y-%m-%d %H:%M:%S', tz='Etc/GMT')

gauge_ids = c(52402.4, 52403.4)
gauge_data = c(
    '../../../../../DATA/TIDES/DART/dart_extract/solomons_2007_04_01_Mw8.1/solomons_2007_04_01_Mw8.1_52402.csv',
    '../../../../../DATA/TIDES/DART/dart_extract/solomons_2007_04_01_Mw8.1/solomons_2007_04_01_Mw8.1_52403.csv')

plot_durations = list(c(0, 3600 * 8), c(0, 3600 * 8))
gauge_ylims = list(c(-1, 1)*0.02, c(-1, 1)*0.02)
event_basename = 'solomon_2007_04_01_Mw81'
## END INPUT PAR

source('check_dart_include.R')


###########################################################################
event_magnitude = 7.8
event_hypocentre = c(161.32, -10.68)
event_start = strptime('2016-12-08 17:38:46', format='%Y-%m-%d %H:%M:%S', tz='Etc/GMT')
gauge_ids = c(55012.4, 55023.4, 52406.4)
gauge_data = c(
    '../../../../../DATA/TIDES/DART/dart_extract/solomons_2016_12_08_Mw7.8/solomons_2016_12_08_Mw7.8_55012.csv',
    '../../../../../DATA/TIDES/DART/dart_extract/solomons_2016_12_08_Mw7.8/solomons_2016_12_08_Mw7.8_55023.csv',
    '../../../../../DATA/TIDES/DART/dart_extract/solomons_2016_12_08_Mw7.8/solomons_2016_12_08_Mw7.8_52406.csv')

plot_durations = list(c(0, 3600 * 6), c(0, 3600 * 6))
gauge_ylims = list(c(-1, 1)*0.05, c(-1, 1)*0.05)
event_basename = 'solomon_2016_12_08_Mw78'
## END INPUT PAR

source('check_dart_include.R')


