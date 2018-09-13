check_dart_env = new.env()
source('compare_with_data_environment.R', local=check_dart_env)

## INPUT PAR
event_magnitude = 7.8
event_hypocentre = c(166.56, -45.76)
event_start = strptime('2009-07-15 09:22:29', format='%Y-%m-%d %H:%M:%S', tz='Etc/GMT')

gauge_ids = c(55015.4, 55042.4)
gauge_data = c(
    '../../../../../DATA/TIDES/DART/dart_extract/puysegur_2009_07_15_Mw7.8/puysegur_2009_07_15_Mw7.8_55015.csv',
    '../../../../../DATA/TIDES/DART/dart_extract/puysegur_2009_07_15_Mw7.8/puysegur_2009_07_15_Mw7.8_55013.csv')

plot_durations = list(c(0, 3600 * 2), c(0, 3600 * 2))
gauge_ylims = list(c(-1, 1)*0.1, c(-1, 1)*0.1)

event_basename = 'puysegur_2009_07_15_Mw7.8'

## END INPUT PAR

source('check_dart_include.R')


