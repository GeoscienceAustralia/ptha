# Key libraries / routines we need
library(rptha)
source('sum_tsunami_unit_sources.R', local=TRUE)
tsunami_model_config = new.env()
source('../TSUNAMI_UNIT_SOURCE/config.R', local=tsunami_model_config, 
    chdir=TRUE)

#
# Data for this source zone
#
source_name = basename(dirname(getwd()))
# NetCDF file with uniform slip earthquake events
earthquake_events = read_table_from_netcdf(
    paste0('all_uniform_slip_earthquake_events_', source_name, '.nc'))
# NetCDF file with stochastic slip earthquake events
earthquake_events_stochastic = read_table_from_netcdf(
    paste0('all_stochastic_slip_earthquake_events_', source_name, '.nc'))
# NetCDF file with stochastic slip earthquake events
earthquake_events_variable_uniform = read_table_from_netcdf(
    paste0('all_variable_uniform_slip_earthquake_events_', source_name, '.nc'))
# NetCDF file with unit-source statistics
unit_source_statistics = read_table_from_netcdf(
    paste0('unit_source_statistics_', source_name, '.nc'))
# Shapefile with unit-source grid geometry
unit_source_geometry = readOGR(dsn='../EQ_SOURCE/unit_source_grid', 
    layer=source_name)
# Tide gauge NetCDF files, sorted based on unit_source_statistics rows
all_tide_files = unit_source_statistics$tide_gauge_file
# Data.frame with the hazard-point (i.e. tide gauge) lon/lat/depth/ID 
all_gauge_lonlat = get_netcdf_gauge_locations(all_tide_files[1]) 
# Get the gauge output times [identical for all files]
gauge_times = get_netcdf_gauge_output_times(all_tide_files[1])


#' Plot uniform slip models and observations
#'
#' @param event_magnitude numeric earthquake magnitude. Events with this
#' magnitude will be used for the plot
#' @param event_hypocentre vector c(lon,lat) giving the location of a point on
#' the rupture. All modelled uniform slip earthquakes will contain a unit-source
#' within 0.5 scaling-law width/length of the unit-source containing this point.
#' (The point must be inside a unit source). If use_stochastic_slip=TRUE, then
#' we will use stochastic events 'corresponding' to the above identified uniform events.
#' @param event_start A POSIX.lt object giving the event start time (UTC), made
#' with e.g.: 
#'     event_start=strptime('2009-06-13 15:22:31', 
#'         format='%Y-%m-%d %H:%M:%S', tz='Etc/UTC')
#' @param gauge_ids vector giving IDs of gauges at which to extract model
#' predictions
#' @param gauge_data character vector of gauge data, with one entry for each
#' gauge_ids. It should give a filename where the de-tided data for that gauge
#' is stored, and be sorted in an order corresponding to gauge_ids.
#' @param plot_durations list with one entry for each gauge_id, containing
#' numeric limits for plot x-axis of the form c(start_sec, end_sec). This gives
#' the plot x-limits in seconds relative to the event_start time. 
#' @param gauge_ylims list with one entry for each gauge_id, containing numeric
#' limits for plot y-axis of the form c(ylower, yupper). This gives the plot
#' y-limits in meters for the de-tided data and the model.
#' @param output_dir_tag character giving a string to insert in the output
#' directory containing model/obs data. If NULL, we do not make those outputs
#' @param use_stochastic_slip logical. If TRUE, generate stochastic slip events.
#' Otherwise use uniform slip
#'
compare_event_with_gauge_time_series<-function(
    event_magnitude, 
    event_hypocentre, 
    event_start, 
    gauge_ids, 
    gauge_data, 
    plot_durations, 
    gauge_ylims, 
    output_dir_tag=NULL,
    use_stochastic_slip = FALSE,
    use_variable_uniform_slip = FALSE){
  
    # Events with the right magnitude. Allow for some rounding error, due to
    # the use of float's in the netcdf file. 
    # Use uniform slip events -- later, if use_stochastic_slip = TRUE, then
    # we will extract those events which correspond to the uniform ones to keep
    keep = which(abs(earthquake_events$Mw - event_magnitude) < 1.0e-03)
    events_with_Mw = earthquake_events[keep, ]
  
    # Find which events contain the hypocentre, by finding which unit source
    # contains it
    unit_source_containing_hypocentre = find_unit_source_index_containing_point(
        event_hypocentre, unit_source_geometry, unit_source_statistics) 

    # Allow events which 'touch' sites within uniform slip scaling law width and half length
    expand_unit_source_alongstrike = ceiling(Mw_2_rupture_size(event_magnitude)[3] * 0.5/
        unit_source_statistics$length[unit_source_containing_hypocentre])
    expand_unit_source_downdip = ceiling(Mw_2_rupture_size(event_magnitude)[2] * 0.5/
        unit_source_statistics$width[unit_source_containing_hypocentre])

    # Find unit-source neighbours (within a few unit-sources in each direction)
    hypocentre_neighbours = c()
    for(j in (seq(-expand_unit_source_downdip,expand_unit_source_downdip))){
        for(i in (seq(-expand_unit_source_alongstrike,expand_unit_source_alongstrike))){
            usch = unit_source_containing_hypocentre # shorthand
            nbr = which(
                (unit_source_statistics$downdip_number == (unit_source_statistics$downdip_number[usch] + j)) &
                (unit_source_statistics$alongstrike_number == (unit_source_statistics$alongstrike_number[usch] + i))
                )
            if(length(nbr) > 1) stop('BUG! This should be impossible')
            if(length(nbr) == 1){
                hypocentre_neighbours = c(hypocentre_neighbours, nbr)
            }
        }
    }
    if(length(hypocentre_neighbours) == 0) stop('No unit source neighbours found')

    event_contains_hpc = rep(0, length(events_with_Mw[,1]))
    for(i in 1:length(events_with_Mw[,1])){
        event_contains_hpc[i] = any(get_unit_source_indices_in_event(
            events_with_Mw[i,]) %in% c(unit_source_containing_hypocentre, hypocentre_neighbours))
    }

    if(sum(event_contains_hpc) == 0){
        stop('No unit sources contain the provided hypocentre')
    }

    if(use_stochastic_slip | use_variable_uniform_slip){

        if(use_variable_uniform_slip & use_stochastic_slip){
            stop('Cannot have both use_stochastic_slip and use_variable_uniform_slip')
        }

        # Find stochastic slip events that correspond to the uniform events we would keep
        # This way of selecting stochastic events should give an unbiased
        # representation of the stochastic model.
        # OTOH, if we just selected 'all stochastic events close enough to the
        # target location', we might introduce bias, since that approach would
        # tend to select more 'broad' events at more distant locations

        if(use_stochastic_slip){
            stochastic_events_uniform_row = 
                earthquake_events_stochastic$uniform_event_row

            uniform_keepers = keep[which(event_contains_hpc == 1)]

            stochastic_events_to_keep = which(
                stochastic_events_uniform_row %in% uniform_keepers)

            # Replace events_with_Mw with the 'stochastic slip' version,
            # containing events that correspond to the uniform events we would
            # keep in the alternate case where use_stochastic_slip=FALSE
            events_with_Mw = earthquake_events_stochastic[stochastic_events_to_keep,]

        }else if(use_variable_uniform_slip){

            variable_uniform_events_uniform_row = 
                earthquake_events_variable_uniform$uniform_event_row

            uniform_keepers = keep[which(event_contains_hpc == 1)]

            variable_uniform_events_to_keep = which(
                variable_uniform_events_uniform_row %in% uniform_keepers)
            # Replace events_with_Mw with the 'variable_uniform slip' version,
            # containing events that correspond to the uniform events we would
            # keep in the alternate case where use_variable_uniform_slip=FALSE
            events_with_Mw = earthquake_events_variable_uniform[
                variable_uniform_events_to_keep,]
        }

    }else{
        # Uniform slip (use_stochastic_slip=FALSE)
        events_with_Mw = events_with_Mw[which(event_contains_hpc == 1),]
    }

    plot_events_vs_gauges(events_with_Mw, event_start, gauge_ids, gauge_data, 
        plot_durations, gauge_ylims, output_dir_tag)

}

#' Stochastic slip variant of \code{compare_event_with_gauge_time_series}
#'
#' @param event_magnitude numeric earthquake magnitude. Make stochastic slip
#' events with this magnitude
#' @param event_hypocentre A vector c(lon,lat) giving the location of a point on
#' the rupture. If create_NEW = true, then all modelled stochastic earthquakes
#' are 'near' this point (peak slip within half-a-width and half-a-length of the
#' location, with width/length based on Strasser earthquake size scaling
#' relations). If create_NEW = FALSE, then the same approach as
#' \code{compare_event_with_gauge_time_series} is used (full scaling law width/length)
#' @param number_of_sffm How many stochastic scenarios to simulate. Beware this
#' is ignored if create_new = FALSE (default)
#' @param zero_low_slip_cells_fraction number close to zero in [0,1). To reduce
#' the number of unit-sources involved in events (and thus reduce the
#' computational and memory requirements), we set slip to zero on cells with
#' smallest slip, which contributing < this fraction of the total cumulative
#' slip. Suggested values of e.g. 0.02 or 0.03. Set to zero to not simplify the
#' slip at all. Beware this is ignored if create_new = FALSE (default)
#' @param event_start POSIX.lt object giving the event start time (UTC), made
#' with e.g. :
#'     event_start=strptime('2009-06-13 15:22:31', format='%Y:%m:%d %H:%M:%S', 
#'         tz='Etc/UTC')
#' @param gauge_ids vector giving IDs of gauges at which to extract model
#' predictions
#' @param gauge_data character vector of gauge data, with one entry for each
#' gauge_ids. It should give a filename where the de-tided data for that gauge
#' is stored, and be ordered to correspond to gauge_ids.
#' @param plot_durations list with one entry for each gauge_id, containing
#' numeric limits for plot x-axis of the form c(start_sec, end_sec). This gives
#' the plot x-limits in seconds relative to the event_start time. 
#' @param gauge_ylims list with one entry for each gauge_id, containing numeric
#' limits for plot y-axis of the form c(ylower, yupper). This gives the plot
#' y-limits in m for the de-tided data and the model.
#' @param output_dir_tag character giving a string to insert in the output
#' directory containing model/obs data. If NULL, we do not make those outputs
#' @param create_new logical. If TRUE, then make new random events. Otherwise,
#' selected events from the existing earthquake_events_stochastic data.frame
#'
compare_stochastic_slip_event_with_gauge_time_series<-function(
    event_magnitude, 
    event_hypocentre, 
    number_of_sffm, 
    zero_low_slip_cells_fraction, 
    event_start, 
    gauge_ids, 
    gauge_data, 
    plot_durations, 
    gauge_ylims,
    output_dir_tag=NULL,
    create_new = FALSE){

    if(create_new){
        # Make new events
        all_events = sffm_make_events_on_discretized_source(
            unit_source_statistics,    
            target_location = event_hypocentre,
            target_event_mw = event_magnitude,
            num_events = number_of_sffm,
            zero_low_slip_cells_fraction=zero_low_slip_cells_fraction,
            sourcename = source_name)

        events_with_Mw = sffm_events_to_table(all_events, 
            slip_significant_figures=4)

        plot_events_vs_gauges(events_with_Mw, event_start, gauge_ids, 
            gauge_data, plot_durations, gauge_ylims, output_dir_tag)
    }else{
        # Get events from existing table
        compare_event_with_gauge_time_series(event_magnitude, event_hypocentre, 
            event_start, gauge_ids, gauge_data, plot_durations, gauge_ylims, 
            output_dir_tag=output_dir_tag,
            use_stochastic_slip = TRUE)
    }

}

#' Plotting code
plot_events_vs_gauges<-function(events_with_Mw, event_start, gauge_ids, 
    gauge_data, plot_durations, gauge_ylims, output_dir_tag=NULL){

    gauge_subset_indices = sapply(gauge_ids, 
        f<-function(x) which.min(abs(x - all_gauge_lonlat$gaugeID)))

    save_outputs = !is.null(output_dir_tag)
    if(save_outputs){
        if('slip' %in% names(events_with_Mw)){
            event_type = 'uniform'
        }else{
            # Slip is either stochastic, or variable_uniform
        
            ess = events_with_Mw$event_slip_string
            ess_split = sapply(ess, f<-function(x) strsplit(x, '_'))
   
            # If all slip values in a single event are the same, we are doing
            # variable uniform 
            is_variable_uniform = 
                all(unlist(lapply(ess_split, f<-function(x) all(x == x[1]))))
    
            if(is_uniform_variable){
                event_type = 'variable_uniform' 
            }else{
                event_type = 'stochastic'
            }
        }
        output_dir = paste0(events_with_Mw$sourcename[1], '_', output_dir_tag, 
            '_', event_type)
        dir.create(output_dir, showWarnings=FALSE)
    }

    for(i in 1:length(gauge_subset_indices)){
        # NOTE: This gets the gauges one at a time. While that's convenient
        # here, it might be faster to do it in blocks ?
        model_events = make_tsunami_event_from_unit_sources(
            events_with_Mw,
            unit_source_statistics,
            all_tide_files,
            indices_of_subset = gauge_subset_indices[i]
            )

        gauge_obs = read.csv(gauge_data[i], stringsAsFactors=FALSE)
        gauge_obs_times = difftime( 
            strptime(gauge_obs$time, format='%Y-%m-%d %H:%M:%S', tz='Etc/GMT'), 
            event_start, 
            units='secs')

        for(j in 1:length(model_events)){
            if(any(is.na(gauge_ylims[[i]]))){
                gauge_ylims[[i]] = c(min(gauge_obs$resid), max(gauge_obs$resid))
            }
            plot(gauge_times, model_events[[j]][1, , 1], t='l', 
                xlab='Time from event (s)', ylab='Stage (m)', 
                xlim=plot_durations[[i]], ylim=c(gauge_ylims[[i]]))
            points(gauge_obs_times, gauge_obs$resid, t='l', col='blue')
            grid()
            abline(v=(-200:200)*3600, col='green')
            title(paste0('Modelled and observed stage, gauge: ', gauge_ids[i], 
                ', Mw: ', events_with_Mw$Mw[j]),
                sub = events_with_Mw$event_index_string[j])
        }

        # Store the model/gauge time-series for other analysis
        if(save_outputs){
            saveRDS(list(model_events = model_events, model_times = gauge_times, 
                gauge_obs = gauge_obs, gauge_obs_times = gauge_obs_times),
                file = paste0(output_dir, '/gauge_', gauge_ids[i], '_Mw_', 
                    events_with_Mw$Mw[1], '.RDS'))
        }

        rm(model_events, gauge_obs, gauge_obs_times)
        gc(verbose=FALSE)
    }

    if(save_outputs){
        saveRDS(list(events_with_Mw = events_with_Mw, 
                unit_source_statistics=unit_source_statistics),
            file=paste0(output_dir, '/event_metadata.RDS'))
    }
}

