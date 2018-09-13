check_dart_env = new.env()
source('compare_with_data_environment.R', local=check_dart_env)

## INPUT PAR
event_magnitude = 7.6
event_hypocentre = c(183.66, -29.54)
event_start = strptime('2011-07-06 19:03:18', format='%Y-%m-%d %H:%M:%S', tz='Etc/GMT')

gauge_ids = c(54401.4)
gauge_data = c(
    '../../../../../DATA/TIDES/DART/dart_extract/kermadectonga_2011_07_06_Mw7.6/kermadectonga_2011_07_06_Mw7.6_54401.csv'
    )

plot_durations = list(c(0, 3600 * 2), c(0, 3600 * 2))
gauge_ylims = list(c(-1, 1)*0.1, c(-1, 1)*0.1)

event_basename = 'outerrise_kermadectonga_2011_07_06_Mw76'

## END INPUT PAR

source('check_dart_include.R')
