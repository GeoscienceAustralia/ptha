library(rptha)
source('sum_tsunami_unit_sources.R', local=TRUE)
tsunami_model_config = new.env()
source('../TSUNAMI_UNIT_SOURCE/config.R', local=tsunami_model_config, chdir=TRUE)

# Data for this source zone
source_name = basename(dirname(getwd()))
#earthquake_events = read.csv(paste0('all_uniform_slip_earthquake_events_', source_name, '.csv'), stringsAsFactors=FALSE)
earthquake_events = read_table_from_netcdf(paste0('all_uniform_slip_earthquake_events_', source_name, '.nc'))
#unit_source_statistics = read.csv(paste0('unit_source_statistics_', source_name, '.csv'), stringsAsFactors=FALSE)
unit_source_statistics = read_table_from_netcdf(paste0('unit_source_statistics_', source_name, '.nc'))
unit_source_geometry = readOGR(dsn='../EQ_SOURCE/unit_source_grid', layer=source_name)
# Get the tide files, sorted based on unit_source_statistics rows
all_tide_files = unit_source_statistics$tide_gauge_file
# Get the gauge lon/lat/depth/ID points
all_gauge_lonlat = get_netcdf_gauge_locations(all_tide_files[1]) 
# Get the gauge output times [identical for all files]
gauge_times = get_netcdf_gauge_output_times(all_tide_files[1])


#' Plot uniform slip models and observations
#'
#' @param event_magnitude numeric earthquake magnitude. Events with this magnitude will
#' be used for the plot
#' @param event_hypocentre. c(lon,lat) giving the location of a point on the
#' rupture. All modelled earthquakes which are plotted will contain this point.
#' The point must be inside a unit source
#' @param event_start POSIX.lt object giving the event start time (UTC), made with e.g. 
#' event_start=strptime('2009-06-13 15:22:31', format='%Y:%m:%d %H:%M:%S', tz='Etc/UTC')
#' @param gauge_ids vector giving IDs of gauges at which to extract model
#' predictions
#' @param gauge_data character vector of gauge data, with one entry for each
#' gauge_ids. It should give a filename where the de-tided data for that gauge
#' is stored. 
#' @param plot_durations list with one entry for each gauge_id, containing numeric limits
#' for plot x-axis of the form c(start_sec, end_sec). This gives the plot x-limits in seconds
#' relative to the event_start time. 
#' @param gauge_ylims list with one entry for each gauge_id, containing numeric limits for
#' plot y-axis of the form c(ylower, yupper). This gives the plot y-limits in m for the de-tided
#' data and the model.
#' @param output_dir_tag character giving a string to insert in the output directory containing
#' model/obs data. If NULL, we do not make those outputs
#'
compare_event_with_gauge_time_series<-function(event_magnitude, event_hypocentre, 
    event_start, gauge_ids, gauge_data, plot_durations, gauge_ylims, output_dir_tag=NULL){
  
    # Events with the right magnitude 
    events_with_Mw = earthquake_events[which(abs(earthquake_events$Mw - event_magnitude) < 1.0e-03), ]
  
    # Find which events contain the hypocentre, by finding which unit source
    # contains it
    unit_source_containing_hypocentre = find_unit_source_index_containing_point(
        event_hypocentre, unit_source_geometry, unit_source_statistics) 

    event_contains_hpc = rep(0, length(events_with_Mw[,1]))
    for(i in 1:length(events_with_Mw[,1])){
        event_contains_hpc[i] = any(get_unit_source_indices_in_event(
            events_with_Mw[i,]) %in% unit_source_containing_hypocentre)
    }

    if(sum(event_contains_hpc) == 0) stop('No unit sources contain the provided hypocentre')

    events_with_Mw = events_with_Mw[which(event_contains_hpc == 1),]

    plot_events_vs_gauges(events_with_Mw, event_start, gauge_ids, gauge_data, 
        plot_durations, gauge_ylims, output_dir_tag)

}

#'
#' Stochastic slip variant of 'compare_event_with_gauge_time_series'
#'
#' @param event_magnitude numeric earthquake magnitude. Make stochastic slip
#' events with this magnitude
#' @param event_hypocentre. c(lon,lat) giving the location of a point on the
#' rupture. All modelled stochastic earthquakes are 'near' this point (peak
#' slip within half-a-width and half-a-length of the location, with width/length
#' based on Strasser earthquake size scaling relations). 
#' @param number_of_sffm How many stochastic scenarios to simulate
#' @param zero_low_slip_cells_fraction number close to zero in [0,1). To reduce the
#' number of unit-sources involved in events (and thus reduce the computational and memory
#' requirements), we set slip to zero on cells with smallest slip, which contributing < this
#' fraction of the total cumulative slip. Suggested values of e.g. 0.02 or 0.03. Set to zero
#' to not simplify the slip at all.
#' @param event_start POSIX.lt object giving the event start time (UTC), made with e.g. 
#' event_start=strptime('2009-06-13 15:22:31', format='%Y:%m:%d %H:%M:%S', tz='Etc/UTC')
#' @param gauge_ids vector giving IDs of gauges at which to extract model
#' predictions
#' @param gauge_data character vector of gauge data, with one entry for each
#' gauge_ids. It should give a filename where the de-tided data for that gauge
#' is stored. 
#' @param plot_durations list with one entry for each gauge_id, containing numeric limits
#' for plot x-axis of the form c(start_sec, end_sec). This gives the plot x-limits in seconds
#' relative to the event_start time. 
#' @param gauge_ylims list with one entry for each gauge_id, containing numeric limits for
#' plot y-axis of the form c(ylower, yupper). This gives the plot y-limits in m for the de-tided
#' data and the model.
#' @param output_dir_tag character giving a string to insert in the output directory containing
#' model/obs data. If NULL, we do not make those outputs
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
    output_dir_tag=NULL){

    all_events = sffm_make_events_on_discretized_source(
        unit_source_statistics,    
        target_location = event_hypocentre,
        target_event_mw = event_magnitude,
        num_events = number_of_sffm,
        zero_low_slip_cells_fraction=zero_low_slip_cells_fraction,
        sourcename = source_name)

    events_with_Mw = sffm_events_to_table(all_events, slip_significant_figures=4)

    plot_events_vs_gauges(events_with_Mw, event_start, gauge_ids, gauge_data, 
        plot_durations, gauge_ylims, output_dir_tag)

}

#' Plotting code
plot_events_vs_gauges<-function(events_with_Mw, event_start, gauge_ids, gauge_data, 
    plot_durations, gauge_ylims, output_dir_tag=NULL){

    gauge_subset_indices = sapply(gauge_ids, f<-function(x) which.min(abs(x - all_gauge_lonlat$gaugeID)))

    save_outputs = !is.null(output_dir_tag)
    if(save_outputs){
        if('slip' %in% names(events_with_Mw)){
            event_type = 'uniform'
        }else{
            event_type = 'stochastic'
        }
        output_dir = paste0(events_with_Mw$sourcename[1], '_', output_dir_tag, '_', event_type)
        dir.create(output_dir, showWarnings=FALSE)
    }

    for(i in 1:length(gauge_subset_indices)){
        # NOTE: This gets the gauges one at a time. While that's convenient here, it might
        # be faster to do it in blocks
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
            plot(gauge_times, model_events[[j]][1, , 1], t='l', xlab='Time from event (s)', 
                ylab='Stage (m)', xlim=plot_durations[[i]], ylim=c(gauge_ylims[[i]]))
            points(gauge_obs_times, gauge_obs$resid, t='l', col='blue')
            grid()
            abline(v=(-200:200)*3600, col='green')
            title(paste0('Modelled and observed stage, gauge: ', gauge_ids[i], ', Mw: ', events_with_Mw$Mw[j]),
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
	#print(memory.profile())
    }

    if(save_outputs){
        saveRDS(list(events_with_Mw = events_with_Mw, unit_source_statistics=unit_source_statistics),
            file=paste0(output_dir, '/event_metadata.RDS'))
    }
}

