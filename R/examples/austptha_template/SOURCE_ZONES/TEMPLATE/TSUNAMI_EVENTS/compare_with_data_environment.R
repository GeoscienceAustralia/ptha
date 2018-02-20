# Key libraries / routines we need
library(rptha)
source('sum_tsunami_unit_sources.R', local=TRUE)
tsunami_model_config = new.env()
source('../TSUNAMI_UNIT_SOURCE/config.R', local=tsunami_model_config, 
    chdir=TRUE)

# Get the NGDC data
ngdc_dir = '../../../../../DATA/TSUNAMI_OBS/NGDC_DATABASE/'
ngdc = new.env()
source(paste0(ngdc_dir, 'ngdc_env.R'), local=ngdc, chdir=TRUE)

# Local config
config_env = new.env()
source('config.R', local=config_env)

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
# NetCDF file with variable-uniform slip earthquake events
earthquake_events_variable_uniform = try(read_table_from_netcdf(
    paste0('all_variable_uniform_slip_earthquake_events_', source_name, '.nc')))
# NetCDF file with unit-source statistics
unit_source_statistics = read_table_from_netcdf(
    paste0('unit_source_statistics_', source_name, '.nc'))
# Tide gauge NetCDF files, sorted based on unit_source_statistics rows
all_tide_files = unit_source_statistics$tide_gauge_file
# Data.frame with the hazard-point (i.e. tide gauge) lon/lat/elev/ID 
all_gauge_lonlat = get_netcdf_gauge_locations(all_tide_files[1]) 
# Get the gauge output times [identical for all files]
gauge_times = get_netcdf_gauge_output_times(all_tide_files[1])
# Shapefile with unit-source grid geometry
unit_source_geometry = readOGR(
    dsn=paste0(dirname(dirname(dirname(unit_source_statistics$initial_condition_file[1]))),
        '/unit_source_grid'),
    layer=source_name)


#' Given hypocentre c(lon,lat), find the index of the unit-source closest to it,
#' and the indices of unit-sources that are within L/2, W/2 of it. Here L,W
#' are computed from the 'median' value of the provided scaling relation
#'
find_unit_sources_near_hypocentre<-function(
    event_hypocentre,
    unit_source_geometry,
    unit_source_statistics,
    event_magnitude){

    # Find which events contain the hypocentre, by finding which unit source
    # contains it
    unit_source_containing_hypocentre = find_unit_source_index_containing_point(
        event_hypocentre, unit_source_geometry, unit_source_statistics) 

    # Allow events which 'touch' sites within uniform slip scaling law width
    # and half length
    local_AWL = Mw_2_rupture_size(event_magnitude, 
        relation=config_env$scaling_relation_type)
    expand_unit_source_alongstrike = ceiling(local_AWL[3] * 0.5/
        unit_source_statistics$length[unit_source_containing_hypocentre])
    expand_unit_source_downdip = ceiling(local_AWL[2] * 0.5/
        unit_source_statistics$width[unit_source_containing_hypocentre])

    # Find unit-source neighbours (within a few unit-sources in each direction)
    hypocentre_neighbours = c()
    for(j in (seq(-expand_unit_source_downdip,expand_unit_source_downdip))){
        for(i in (seq(-expand_unit_source_alongstrike,expand_unit_source_alongstrike))){
            usch = unit_source_containing_hypocentre # shorthand
            nbr = which(
                (unit_source_statistics$downdip_number == 
                    (unit_source_statistics$downdip_number[usch] + j)) &
                (unit_source_statistics$alongstrike_number == 
                    (unit_source_statistics$alongstrike_number[usch] + i))
                )
            if(length(nbr) > 1) stop('BUG! This should be impossible')
            if(length(nbr) == 1){
                hypocentre_neighbours = c(hypocentre_neighbours, nbr)
            }
        }
    }
    if(length(hypocentre_neighbours) == 0) stop('No unit source neighbours found')

    output = c(unit_source_containing_hypocentre, hypocentre_neighbours)

    return(output)

}

#'
#' Find earthquake events near a point close to a given magnitude
#'
#' This is actually just a wrapper for other functions which treat
#' the 'variable mu' and 'fixed mu' cases. See those functions for
#' documentation
#'  
#'
find_events_near_point<-function(
    event_magnitude,
    event_hypocentre,
    use_stochastic_slip_runs = FALSE,
    use_variable_uniform_slip_runs = FALSE,
    fixed_mu = TRUE){

    # Local input args
    event_magnitude1 = event_magnitude
    event_hypocentre1 = event_hypocentre
    use_stochastic_slip_runs1 = use_stochastic_slip_runs
    use_variable_uniform_slip_runs1 = use_variable_uniform_slip_runs
    fixed_mu1 = fixed_mu

    if(fixed_mu1){
        # Use fixed mu. 
        events_with_Mw = find_events_near_point_fixed_mu(
            event_magnitude = event_magnitude1,
            event_hypocentre = event_hypocentre1,
            use_stochastic_slip_runs = use_stochastic_slip_runs1,
            use_variable_uniform_slip_runs = use_variable_uniform_slip_runs1)
    }else{
        # Use variable mu
        events_with_Mw = find_events_near_point_variable_mu(
            event_magnitude = event_magnitude1,
            event_hypocentre = event_hypocentre1,
            use_stochastic_slip_runs = use_stochastic_slip_runs1,
            use_variable_uniform_slip_runs = use_variable_uniform_slip_runs1)
    }

    return(events_with_Mw)
}


#' Function similar to find_events_near_point_fixed_mu, but with variable mu.
#' i.e. considering shear modulus depth dependence, re-compute the magnitude
#' of every event and find events 'near' the desired magnitude that have
#' unit-sources near the desired hypocentre
#'
#' @param event_magnitude numeric earthquake magnitude. Events close to this
#' magnitude will be selected 
#' @param event_hypocentre vector c(lon,lat) giving the location of a point on
#' the rupture. All modelled uniform slip earthquakes will contain a unit-source
#' within 0.5 scaling-law width/length of the unit-source containing this point.
#' (The point must be inside a unit source). If use_stochastic_slip_runs=TRUE,
#' then we will use stochastic events 'corresponding' to the above identified
#' uniform events.  If use_variable_uniform_slip_runs = TRUE, then we do as for
#' stochastic_slip, in the variable_uniform case.
#' @param use_stochastic_slip_runs logical. If TRUE, return stochastic slip
#' events that 'correspond to' the uniform slip event
#' @param use_variable_uniform_slip_runs logical. If TRUE, return variable
#' uniform slip events that 'correspond to' the uniform slip event
#' @return data.frame with metadata for the relevant events
#'
find_events_near_point_variable_mu<-function(
    event_magnitude,
    event_hypocentre,
    use_stochastic_slip_runs = FALSE,
    use_variable_uniform_slip_runs = FALSE){

    if(use_stochastic_slip_runs & use_variable_uniform_slip_runs){
        stop("Cannot use both stochastic and variable slip at once")
    }

    # Compute shear modulus given depth
    mu_fun<-function(dpth){
        ## Point values based on plot of Lay and Bilek (2007), limited to 10GPA in shallow areas
        ## See
        #/media/gareth/Windows7_OS/Users/gareth/Documents/work/AustPTHA/DATA/EARTHQUAKE/Shear_modulus
        depths = c(0, 7.5, 15, 35, 9999)
        mu = c(10, 10, 30, 67, 67)*1e+09
        output = 10**(approx(depths, log10(mu), xout=dpth)$y)
        return(output)
    }

    #
    # Compute Mw, considering variation in shear modulus
    #

    # Preliminary variables from unit-source statistics
    depth = unit_source_statistics$depth
    mu_variable = mu_fun(depth)
    area = unit_source_statistics$length*unit_source_statistics$width

    # Get event indices, event slips, and mus
    if(use_stochastic_slip_runs | use_variable_uniform_slip_runs){
        # Not uniform slip format

        if(use_stochastic_slip_runs){
            events = earthquake_events_stochastic
        }else if(use_variable_uniform_slip_runs){
            events = earthquake_events_variable_uniform
        }          

        event_inds = sapply(events$event_index_string, 
            f<-function(x) as.numeric(strsplit(x, '-')[[1]]), simplify=FALSE)
        event_slips = sapply(events$event_slip_string, 
            f<-function(x) as.numeric(strsplit(x, '_')[[1]]), simplify=FALSE)

        # Get the indices in the 'corresponding uniform slip' event
        event_inds_uniform = sapply(earthquake_events$event_index_string,
            f<-function(x) as.numeric(strsplit(x, '-')[[1]]), simplify=FALSE)
        event_inds_uniform = event_inds_uniform[events$uniform_event_row]

    }else{
        # Uniform slip format. We extract data in the same form as would be used for variable slip
        events = earthquake_events

        event_inds = sapply(events$event_index_string, 
            f<-function(x) as.numeric(strsplit(x, '-')[[1]]), simplify=FALSE)
        # Expand slip to the same format as we have for variable slip events
        event_slips = event_inds
        for(i in 1:length(event_slips)){
            event_slips[[i]] = event_slips[[i]] * 0 + events$slip[i]
        }

        # Get the indices in the 'corresponding uniform slip' event
        # Obviously for uniform slip this is just 'event_inds', but for
        # other slip styles it is different
        event_inds_uniform = event_inds

    }

    # Seismic moment
    moment = rep(NA, length(event_inds))
    for(i in 1:length(event_inds)){
        slip = event_slips[[i]]
        inds = event_inds[[i]]
        areas = area[inds]
        mu = mu_variable[inds]

        moment[i] = sum(slip * areas * 1e+06 * mu)
    }

    event_Mw_variable_mu = M0_2_Mw(moment)
    # Store magnitude in the outputs
    events$Mw_variable_mu = event_Mw_variable_mu

    # Find indices of unit sources containing the hypocentre, and those within
    # 1/2 length and 1/2 width of an earthquake with typical size (given the scaling 
    # relation)
    unit_sources_near_hypocentre = find_unit_sources_near_hypocentre(
        event_hypocentre,
        unit_source_geometry,
        unit_source_statistics,
        event_magnitude)

    # Find events with the right magnitude, and which involve the right unit sources
    keep = rep(FALSE, length(moment))
    for(i in 1:length(keep)){
        test1 = ( abs(event_Mw_variable_mu[i] - event_magnitude) < (config_env$dMw*1.5) )
        # Keep events if the CORRESPONDING UNIFORM EVENT contained any unit sources near the hypocentre
        # We do this to avoid bias due to the random event size. If we performed a similar test
        # based on 'event_inds', then large, low slip events would have a disproportionately high
        # chance of being selected because they have a larger area. This would distort the comparison
        test2 = any( event_inds_uniform[[i]] %in% unit_sources_near_hypocentre)
        if(test1 & test2) keep[i] = TRUE 
    }

    if(sum(keep) > 0){
        output = events[which(keep),]
    }else{
        stop('No events')
    }

    return(output)

}

#' Given a magnitude and hypocentre, find earthquake events with the 'same'
#' magnitude which contain (or are near) the hypocentre
#'
#' @param event_magnitude numeric earthquake magnitude. Events with this
#' magnitude will be used for the plot
#' @param event_hypocentre vector c(lon,lat) giving the location of a point on
#' the rupture. All modelled uniform slip earthquakes will contain a unit-source
#' within 0.5 scaling-law width/length of the unit-source containing this point.
#' (The point must be inside a unit source). If use_stochastic_slip_runs=TRUE,
#' then we will use stochastic events 'corresponding' to the above identified
#' uniform events.  If use_variable_uniform_slip_runs = TRUE, then we do as for
#' stochastic_slip, in the variable_uniform case.
#' @param use_stochastic_slip_runs logical. If TRUE, return stochastic slip
#' events that 'correspond to' the uniform slip event
#' @param use_variable_uniform_slip_runs logical. If TRUE, return variable
#' uniform slip events that 'correspond to' the uniform slip event
#' @return data.frame with metadata for the relevant events
find_events_near_point_fixed_mu<-function(
    event_magnitude,
    event_hypocentre,
    use_stochastic_slip_runs = FALSE,
    use_variable_uniform_slip_runs = FALSE){

    # Events with the right magnitude. Allow for some rounding error, due to
    # the use of float's in the netcdf file. Also allow to select events that are
    # 'within dMw*1.5' of the desired magnitude, which is relevant for variable mu cases.
    # Use uniform slip events -- later, if use_stochastic_slip = TRUE, then
    # we will extract those events which correspond to the uniform ones to keep
    keep = which(abs(earthquake_events$Mw - event_magnitude) < (config_env$dMw*1.5) )
    events_with_Mw = earthquake_events[keep, ]
 
    # Find indices of unit sources containing the hypocentre, and those within
    # 1/2 length and 1/2 width of an earthquake with typical size given the scaling 
    # relation
    unit_sources_near_hypocentre = find_unit_sources_near_hypocentre(
        event_hypocentre,
        unit_source_geometry,
        unit_source_statistics,
        event_magnitude)
    
    event_contains_hpc = rep(0, length(events_with_Mw[,1]))
    for(i in 1:length(events_with_Mw[,1])){
        event_contains_hpc[i] = any(
            get_unit_source_indices_in_event(events_with_Mw[i,]) %in% 
            unit_sources_near_hypocentre
            )
    }

    if(sum(event_contains_hpc) == 0){
        stop('No unit sources contain the provided hypocentre')
    }

    if(use_stochastic_slip_runs | use_variable_uniform_slip_runs){

        if(use_variable_uniform_slip_runs & use_stochastic_slip_runs){
            stop('Cannot have both use_stochastic_slip_runs and use_variable_uniform_slip_runs')
        }

        # Find stochastic slip events that correspond to the uniform events we
        # would keep
        # This way of selecting stochastic events should give an unbiased
        # representation of the stochastic model.
        # OTOH, if we just selected 'all stochastic events close enough to the
        # target location', we might introduce bias, since that approach would
        # tend to select more 'broad' events at more distant locations

        uniform_keepers = keep[which(event_contains_hpc == 1)]

        if(use_stochastic_slip_runs){

            stochastic_events_uniform_row = 
                earthquake_events_stochastic$uniform_event_row

            stochastic_events_to_keep = which(
                stochastic_events_uniform_row %in% uniform_keepers)

            # Replace events_with_Mw with the 'stochastic slip' version,
            # containing events that correspond to the uniform events we would
            # keep in the alternate case where use_stochastic_slip_runs=FALSE
            events_with_Mw = earthquake_events_stochastic[stochastic_events_to_keep,]

        }else if(use_variable_uniform_slip_runs){

            variable_uniform_events_uniform_row = 
                earthquake_events_variable_uniform$uniform_event_row

            variable_uniform_events_to_keep = which(
                variable_uniform_events_uniform_row %in% uniform_keepers)
            # Replace events_with_Mw with the 'variable_uniform slip' version,
            # containing events that correspond to the uniform events we would
            # keep in the alternate case where
            # use_variable_uniform_slip_runs=FALSE
            events_with_Mw = earthquake_events_variable_uniform[
                variable_uniform_events_to_keep,]
        }

    }else{
        # Uniform slip
        events_with_Mw = events_with_Mw[which(event_contains_hpc == 1),]
    }

    return(events_with_Mw)

}


#' Plot uniform (or stochastic or uniform-variable) slip models and observations
#'
#' @param event_magnitude numeric earthquake magnitude. Events with this
#' magnitude will be used for the plot
#' @param event_hypocentre vector c(lon,lat) giving the location of a point on
#' the rupture. All modelled uniform slip earthquakes will contain a unit-source
#' within 0.5 scaling-law width/length of the unit-source containing this point.
#' (The point must be inside a unit source). If use_stochastic_slip=TRUE, then
#' we will use stochastic events 'corresponding' to the above identified
#' uniform events.
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
#' @param use_stochastic_slip logical. If TRUE, return stochastic slip events
#' that 'correspond to' the uniform slip event
#' @param use_variable_uniform_slip logical. If TRUE, return variable uniform
#' slip events that 'correspond to' the uniform slip event
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
    use_variable_uniform_slip = FALSE,
    make_plot=TRUE,
    fixed_mu=TRUE){

    make_plot = make_plot

    events_with_Mw = find_events_near_point(
        event_magnitude, 
        event_hypocentre,
        use_stochastic_slip_runs = use_stochastic_slip,
        use_variable_uniform_slip_runs = use_variable_uniform_slip,
        fixed_mu = fixed_mu)  

    plot_events_vs_gauges(events_with_Mw, event_start, gauge_ids, gauge_data, 
        plot_durations, gauge_ylims, output_dir_tag, make_plot=make_plot)

}

#' Stochastic slip variant of \code{compare_event_with_gauge_time_series}
#'
#' Note this routine is largely defunct, because \code{compare_event_with_gauge_time_series}
#' treats all cases [so long as tsunami events have already been created]. This
#' routine cannot treat variable 'mu' (because that is a currently post-processing step,
#' not deep in the rptha package)
#'
#' @param event_magnitude numeric earthquake magnitude. Make stochastic slip
#' events with this magnitude
#' @param event_hypocentre A vector c(lon,lat) giving the location of a point on
#' the rupture. If create_NEW = TRUE, then all modelled stochastic earthquakes
#' are 'near' this point (peak slip within half-a-width and half-a-length of the
#' location, with width/length based on Strasser earthquake size scaling
#' relations). If create_NEW = FALSE, then the same approach as
#' \code{compare_event_with_gauge_time_series} is used (finding 'nearby' uniform slip
#' events with fixed dimensions, and keeping corresponding stochastic/variable_uniform
#' events)
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
    create_new = FALSE,
    make_plot = TRUE){

    make_plot = make_plot

    if(create_new){
        # Make new events
        all_events = sffm_make_events_on_discretized_source(
            unit_source_statistics,    
            target_location = event_hypocentre,
            target_event_mw = event_magnitude,
            num_events = number_of_sffm,
            zero_low_slip_cells_fraction=zero_low_slip_cells_fraction,
            sourcename = source_name,
            mu=config_env$shear_modulus,
            relation=config_env$scaling_relation_type)

        events_with_Mw = sffm_events_to_table(all_events, 
            slip_significant_figures=4)

        plot_events_vs_gauges(events_with_Mw, event_start, gauge_ids, 
            gauge_data, plot_durations, gauge_ylims, output_dir_tag,
            make_plot=make_plot)
    }else{
        # Get events from existing table
        compare_event_with_gauge_time_series(event_magnitude, event_hypocentre, 
            event_start, gauge_ids, gauge_data, plot_durations, gauge_ylims, 
            output_dir_tag=output_dir_tag,
            use_stochastic_slip = TRUE,
            make_plot = make_plot)
    }

}

#'
#' Get modelled peak-stage data for comparison with ncdc
#'
#'
compare_event_maxima_with_NGDC<-function(
    start_date, 
    event_Mw, 
    event_hypocentre,
    use_stochastic_slip = FALSE, 
    use_variable_uniform_slip = FALSE,
    output_dir_tag=NULL,
    fixed_mu = TRUE){

    # Get the NGDC data for the event (on the same day)
    event_year = as.numeric(format(start_date, '%Y'))
    event_month = as.numeric(format(start_date, '%m'))
    event_day = as.numeric(format(start_date, '%d'))

    tsunami_obs = ngdc$get_tsunami_data_on_date(event_year,
        event_month, event_day)

    if(nrow(tsunami_obs) == 0){
        print('No tsunami data found in NGDC database')
        return(invisible())
    }


    # Get the database events that are similar
    events_with_Mw = find_events_near_point(
        event_Mw,
        event_hypocentre,
        use_stochastic_slip_runs = use_stochastic_slip,
        use_variable_uniform_slip_runs = use_variable_uniform_slip,
        fixed_mu = fixed_mu
        )

    matching_event_rows = as.numeric(rownames(events_with_Mw))

    #
    # Open netcdf with peak stage info
    #
    
    # Default case (maybe overwritten below)
    peak_stage_ncdf = paste0('all_uniform_slip_earthquake_events_tsunami_', 
        source_name, '.nc')

    # Fix the name, if stochastic or uniform is used
    if(use_stochastic_slip){

        stopifnot(use_variable_uniform_slip == FALSE)

        peak_stage_ncdf = paste0(
            'all_stochastic_slip_earthquake_events_tsunami_', 
            source_name, '.nc')   
    }    
    if(use_variable_uniform_slip){

        stopifnot(use_stochastic_slip == FALSE)

        peak_stage_ncdf = paste0(
            'all_variable_uniform_slip_earthquake_events_tsunami_', 
            source_name, '.nc')   
    }
    fid = nc_open(peak_stage_ncdf, readunlim=FALSE)

    # Find stations matching each observation
    matching_stations = lonlat_nearest_neighbours(
        cbind(tsunami_obs$LONGITUDE, tsunami_obs$LATITUDE), 
        all_gauge_lonlat[,1:2])

    # Get peak stages for all events and all gauges with observations
    max_stage_store = matrix(NA, nrow=length(matching_event_rows), 
        ncol = length(matching_stations))

    for(i in 1:length(matching_stations)){

        all_peak_stage = ncvar_get(fid, 'max_stage', 
            start = c(1, matching_stations[i]), 
            count=c(-1,1))

        max_stage_store[,i] = all_peak_stage[matching_event_rows]

    }

    # Green's law correction
    correction_factors = (-all_gauge_lonlat$elev[matching_stations])**0.25
    # No correction at DART buoys -- assume gauge already has same depth
    # (will be OK so long as gauge is close to DART)
    kk = which(tsunami_obs$TYPE_MEASUREMENT_ID == 3)
    correction_factors[kk] = 1
    # Assume tide gauges in 10m water depth (for reference comparison purposes)
    # Experience indicates that tide gauges are generally lower than '1m-depth
    # greens law'
    kk = which(tsunami_obs$TYPE_MEASUREMENT_ID == 2)
    correction_factors[kk] = 
        (-all_gauge_lonlat$elev[matching_stations[kk]]/10)**0.25

    # Record the distance between observation point and nearest model gauge.
    # It is wise to filter-out pairs that are not sufficiently close.
    distance_to_nearest = distHaversine(
        cbind(tsunami_obs$LONGITUDE, tsunami_obs$LATITUDE), 
        all_gauge_lonlat[matching_stations,1:2])

    output = list(
        events = events_with_Mw,
        matching_stations = matching_stations,
        gauge_lonlat = all_gauge_lonlat[matching_stations,],
        max_stage = max_stage_store,
        greens_correction_factors=correction_factors,
        gcf_mat = matrix(correction_factors, ncol=ncol(max_stage_store), 
            nrow=nrow(max_stage_store), byrow=TRUE),
        distance_to_nearest = distance_to_nearest,
        tsunami_obs = tsunami_obs,
        start_date = start_date,
        event_Mw = event_Mw,
        event_hypocentre = event_hypocentre)

    if(!is.null(output_dir_tag)){
        # Make output directory if it does not already exist
        # NOTE THE DIRECTORY NAME IS THE SAME AS IN \code{plot_events_vs_gauges}
        # THAT'S IMPORTANT! FIXME: Clean this up!
        event_type = 'uniform'
        if(use_variable_uniform_slip){
            event_type = 'variable_uniform' 
        }
        if(use_stochastic_slip){
            event_type = 'stochastic'
        }
        output_dir = paste0(events_with_Mw$sourcename[1], '_', output_dir_tag, 
            '_', event_type)
        dir.create(output_dir, showWarnings=FALSE)

        saveRDS(output, paste0(output_dir, '/event_NGDC_comparison.RDS'))
    }

    return(invisible(output))

}

#'
#' Plotting code for gauges. Also saves outputs in RDS format (important
#' side-effect that is exploited by other code in ./plots). In some instances
#' we might not want the plot, just the latter side effect, and can pass
#' make_plot=FALSE to ensure that.
#'
plot_events_vs_gauges<-function(events_with_Mw, event_start, gauge_ids, 
    gauge_data, plot_durations, gauge_ylims, output_dir_tag=NULL,
    make_plot=TRUE){

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
    
            if(is_variable_uniform){
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

        if(make_plot){
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
                if('Mw_variable_mu' %in% names(events_with_Mw)){
                    mw_name = round(events_with_Mw$Mw_variable_mu[j], 3)
                }else{
                    mw_name = round(events_with_Mw$Mw[j], 3)
                }
                title(paste0('Modelled and observed stage, gauge: ', gauge_ids[i], 
                    ', Mw: ', mw_name),
                    sub = events_with_Mw$event_index_string[j])
            }
        }

        # Store the model/gauge time-series for other analysis
        if(save_outputs){
            if('Mw_variable_mu' %in% names(events_with_Mw)){
                mw_name = round(median(events_with_Mw$Mw_variable_mu), 3)
            }else{
                mw_name = round(median(events_with_Mw$Mw), 3)
            }

            saveRDS(list(model_events = model_events, model_times = gauge_times, 
                gauge_obs = gauge_obs, gauge_obs_times = gauge_obs_times),
                file = paste0(output_dir, '/gauge_', gauge_ids[i], '_Mw_', 
                    mw_name, '.RDS'))
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

