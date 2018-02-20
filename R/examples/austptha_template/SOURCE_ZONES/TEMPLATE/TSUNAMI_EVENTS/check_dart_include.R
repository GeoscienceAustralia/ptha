#
# This code can be 'included' or 'sourced' in a check_dart.R script that
# defines:
#  event_magnitude, event_hypocentre, gauge_ids, gauge_data, plot_durations, gauge_ylims, event_basename
#
# Useful to reduce code replication
#

stopifnot(exists('event_basename'))
#
# Uniform slip
#

for(variable_mu in c(TRUE, FALSE)){

    if(variable_mu){
        extra_name = '_varyMu'
    }else{
        extra_name = ''
    }

    check_dart_env$compare_event_with_gauge_time_series(
        event_magnitude,
        event_hypocentre,
        event_start,
        gauge_ids,
        gauge_data,
        plot_durations= plot_durations,
        gauge_ylims = gauge_ylims,
        output_dir_tag = paste0(event_basename, extra_name, '_uniform'),
        make_plot=FALSE,
        fixed_mu = (!variable_mu))
    # Store NGDC comparison data
    check_dart_env$compare_event_maxima_with_NGDC(
        event_start,
        event_magnitude,
        event_hypocentre,
        output_dir_tag = paste0(event_basename, extra_name, '_uniform'),
        fixed_mu = (!variable_mu))

    #
    # Stochastic slip
    #

    check_dart_env$compare_event_with_gauge_time_series(
        event_magnitude,
        event_hypocentre,
        event_start,
        gauge_ids,
        gauge_data,
        plot_durations= plot_durations,
        gauge_ylims = gauge_ylims,
        output_dir_tag = paste0(event_basename, extra_name, '_stochastic'),
        make_plot=FALSE,
        use_stochastic_slip=TRUE,
        fixed_mu = (!variable_mu))
    # Store NGDC comparison data
    check_dart_env$compare_event_maxima_with_NGDC(
        event_start,
        event_magnitude,
        event_hypocentre,
        output_dir_tag = paste0(event_basename, extra_name, '_stochastic'),
        use_stochastic_slip = TRUE,
        fixed_mu = (!variable_mu))

    #
    # Variable uniform slip
    #

    check_dart_env$compare_event_with_gauge_time_series(
        event_magnitude,
        event_hypocentre,
        event_start,
        gauge_ids,
        gauge_data,
        plot_durations= plot_durations,
        gauge_ylims = gauge_ylims,
        output_dir_tag = paste0(event_basename, extra_name, '_variable_uniform'),
        make_plot=FALSE,
        use_variable_uniform=TRUE,
        fixed_mu = (!variable_mu))
    # Store NGDC comparison data
    check_dart_env$compare_event_maxima_with_NGDC(
        event_start,
        event_magnitude,
        event_hypocentre,
        output_dir_tag = paste0(event_basename, extra_name, '_variable_uniform'),
        use_variable_uniform_slip = TRUE,
        fixed_mu = (!variable_mu))

}
