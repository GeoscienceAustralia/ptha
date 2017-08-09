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

## END INPUT PAR

#
# Uniform slip
#

pdf('puysegur_2009_07_15_Mw78.pdf', width=10, height=5)
check_dart_env$compare_event_with_gauge_time_series(
    event_magnitude,
    event_hypocentre,
    event_start,
    gauge_ids,
    gauge_data,
    plot_durations= plot_durations,
    gauge_ylims = gauge_ylims,
    output_dir_tag = 'puysegur_2009_07_15_Mw78_uniform')
dev.off()
# Store NGDC comparison data
compare_event_maxima_with_NGDC(
    event_start,
    event_magnitude,
    event_hypocentre,
    output_dir_tag = 'puysegur_2009_07_15_Mw78_uniform')

#
# Stochastic slip
#

pdf('puysegur_2009_07_15_Mw78_stochastic.pdf', width=10, height=5)
check_dart_env$compare_event_with_gauge_time_series(
    event_magnitude,
    event_hypocentre,
    event_start,
    gauge_ids,
    gauge_data,
    plot_durations= plot_durations,
    gauge_ylims = gauge_ylims,
    output_dir_tag = 'puysegur_2009_07_15_Mw78_stochastic',
    use_stochastic_slip=TRUE)
dev.off()
# Store NGDC comparison data
compare_event_maxima_with_NGDC(
    event_start,
    event_magnitude,
    event_hypocentre,
    output_dir_tag = 'puysegur_2009_07_15_Mw78_stochastic',
    use_stochastic_slip = TRUE)

#
# Variable uniform slip
#

pdf('puysegur_2009_07_15_Mw78_variable_uniform.pdf', width=10, height=5)
check_dart_env$compare_event_with_gauge_time_series(
    event_magnitude,
    event_hypocentre,
    event_start,
    gauge_ids,
    gauge_data,
    plot_durations= plot_durations,
    gauge_ylims = gauge_ylims,
    output_dir_tag = 'puysegur_2009_07_15_Mw78_variable_uniform',
    use_variable_uniform=TRUE
    )
dev.off()
# Store NGDC comparison data
compare_event_maxima_with_NGDC(
    event_start,
    event_magnitude,
    event_hypocentre,
    output_dir_tag = 'puysegur_2009_07_15_Mw78_variable_uniform',
    use_variable_uniform_slip = TRUE)


