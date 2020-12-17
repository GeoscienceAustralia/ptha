#
# Codes to remotely access various datasets from new Australian PTHA, on the NCI
#

suppressPackageStartupMessages(library(raster))
if(!exists('config_env')){
    config_env = new.env()
    source('R/config.R', local=config_env)
}
get_supporting_data = new.env()
source('R/get_supporting_data.R', local=get_supporting_data)
source('R/sum_tsunami_unit_sources.R', local=TRUE)
 

#' Read key summary statistics for earthquake events on the source-zone
#'
#' The code reads the files from the web, so an internet connection is required.
#' At the moment no caching of data is implemented. Also, the file sizes range
#' greatly (from a few MB up to 15GB). Subsetting must be used on the large files
#' (assuming typical 2018 internet speeds).
#'
#' @param source_zone Name of source_zone
#' @param slip_type 'stochastic' or 'variable_uniform' or 'uniform'
#' @param desired_event_rows integer vector giving the rows of the table that
#' are desired. If NULL, read all rows (unless range_list is not NULL, see below)
#' @param range_list If not NULL, this list is used for selecting subsets of the
#' events data. For example if range_list=list(Mw=c(9.05, 9.15), peak_slip_alongstrike_ind=c(80,90)),
#' then the selected events will all have Mw in >9.05 and <9.15, and peak_slip_alongstrike_ind >80 and < 90.
#' If range_list is specified, you should not provide desired_event_rows.
#' @param chunk_size The chunk_size passed to read_table_from_netcdf. This can impact the efficiency when only
#' reading a subset of the file. A small chunk size will lead to many reads, whereas a large chunk size may lead
#' to too much data being read at once. The best size depends on the system and the interconnect.
#' @return list with 'events' giving summary statistics for the earthquake
#' events, and 'unit_source_statistics' giving summary statistics for each
#' unit source, and 'gauge_netcdf_files' giving the tide-gauge netcdf filenames for each unit_source,
#' and 'desired_event_rows' giving the row_indices that were selected from the full events table file,
#' and 'slip_type'
#' @export
#' @examples
#' # Basic usage
#' puysegur_data = get_source_zone_events_data('puysegur', slip_type='stochastic')
#'
#' # Select only event_table rows in 10-20 and 40-45
#' puysegur_data_subset = get_source_zone_events_data('puysegur', desired_event_rows = c(10:20, 40:45))
#'
#' # Use a range_list to only select events with Mw~8.1, with rate_annual > 0
#' puysegur_data_Mw81 = get_source_zone_events_data('puysegur', range_list=list(Mw=c(8.05, 8.15), rate_annual=c(0, Inf)))
#'
get_source_zone_events_data<-function(source_zone=NULL, slip_type='stochastic', desired_event_rows = NULL,
                                      range_list=NULL, chunk_size=1000, include_potential_energy=FALSE){

    library(rptha)

    # First check that a valid source-zone was provided
    err = FALSE
    if(is.null(source_zone)){
        err = TRUE
    }else{
        if(sum(config_env$source_names_all == source_zone) == 0) err=TRUE
    }

    if(err){
        print('You did not pass a valid source_zone to get_source_zone_events_data. The allowed source_zone values are:')
        print(paste0('   ', config_env$source_names_all))
        print('Please pass one of the above source_zone names to this function to get its metadata')
        # Fail gracefully
        output = list(events = NA, unit_source_statistics=NA, gauge_netcdf_files=NA)
        return(invisible(output))
    }

    stopifnot(slip_type %in% c('uniform', 'stochastic', 'variable_uniform'))

    #
    # Get the earthquake events data
    #
    nc_web_addr = paste0(config_env$.GDATA_OPENDAP_BASE_LOCATION, 
        'SOURCE_ZONES/', source_zone, '/TSUNAMI_EVENTS/all_', slip_type, 
        '_slip_earthquake_events_', source_zone, '.nc')    

    # Identify desired_event_rows using the range_list, if provided
    if(!is.null(range_list)){
        if(!is.null(desired_event_rows)) stop('Cannot provide BOTH desired_event_rows AND range_list')
        var_data = read_table_from_netcdf(nc_web_addr, varnames=names(range_list))
        to_keep = rep(TRUE, length(var_data[,1])) # Predefine
        for(i in 1:length(range_list)){
            vname = names(range_list)[i]
            to_keep = (to_keep & 
                       (var_data[[vname]] >= range_list[[i]][1]) & 
                       (var_data[[vname]] <= range_list[[i]][2]) )
            if(!any(to_keep)) stop(paste0('No events within the range_list bounds. Failed at variable_name=', vname))
        }
        desired_event_rows = which(to_keep)
    }

    events_file = nc_web_addr
    events_data = read_table_from_netcdf(events_file, desired_rows = desired_event_rows, chunk_size=chunk_size)

    if(include_potential_energy){
        # Read the initial condition potential energy information as well
        potential_energy_data = get_source_zone_events_potential_energy(source_zone, 
            slip_type, desired_event_rows, chunk_size)
        # Logical checks
        if(nrow(potential_energy_data) != nrow(events_data)){
            stop('error reading potential energy -- not aligned with event_data')
        }

        # More sanity checks. 
        # Confirm that Mw, rate_annual, and weight_with_nonzero_rate are the same
        t1 = all(abs(events_data$Mw - potential_energy_data$Mw) <= 1.e-04)
        t2 = all(abs(events_data$rate_annual - potential_energy_data$rate_annual) <= 
                 1.0e-04*events_data$rate_annual)
        t3 = all(abs(events_data$weight_with_nonzero_rate - potential_energy_data$weight_with_nonzero_rate) <=
                 1.0e-04*events_data$weight_with_nonzero_rate)

        if(!(t1 & t2 & t3)){
            stop('Error reading potential energy -- rows may not be aligned')
        }

        # If we got here, then it all worked
        events_data$initial_potential_energy = potential_energy_data$initial_potential_energy
        rm(potential_energy_data)
        gc()
    }

    # Record the tsunami events file too, although we don't use it here
    tsunami_events_file = paste0(config_env$.GDATA_OPENDAP_BASE_LOCATION, 
        'SOURCE_ZONES/', source_zone, '/TSUNAMI_EVENTS/all_', slip_type, 
        '_slip_earthquake_events_tsunami_', source_zone, '.nc')

    #
    # Get the unit source summary statistics
    #
    nc_web_addr = paste0(config_env$.GDATA_OPENDAP_BASE_LOCATION, 'SOURCE_ZONES/', 
        source_zone, '/TSUNAMI_EVENTS/unit_source_statistics_', source_zone, 
        '.nc')    

    unit_source_file = nc_web_addr
    unit_source_statistics = read_table_from_netcdf(unit_source_file)

    gauge_netcdf_files = unit_source_statistics$tide_gauge_file 

    # Append the web address to the files 
    gauge_netcdf_files_reordered = sapply(gauge_netcdf_files, 
        config_env$adjust_path_to_gdata_base_location)

    output = list()
    output[['events']] = events_data
    output[['unit_source_statistics']] = unit_source_statistics
    output[['gauge_netcdf_files']] = gauge_netcdf_files_reordered
    output[['desired_event_rows']] = desired_event_rows
    output[['events_file']] = events_file
    output[['unit_source_file']] = unit_source_file
    output[['tsunami_events_file']] = tsunami_events_file

    return(invisible(output))
}

#' Create initial conditions (i.e. water surface deformation) for tsunami model
#'
#' @param source_zone_events_data output from \code{get_source_zone_events_data}
#' @param event_ID integer ID, corresponding to the row-index of the event in
#' the events table (which is contained in source_zone_events_data) 
#' @param force_file_download logical. If FALSE, we check for files in the local
#' cache and use those if they exist, or download them otherwise. If TRUE, we 
#' download the files to the local cache, irrespective of whether files with
#' the same name already exist.
#' @return A raster with the initial free surface deformation
#' @export
#' @examples
#' puysegur_data = get_source_zone_events_data('puysegur')
#' # Get initial condition corresponding to the event in
#' # puysegur$events[250,]
#' initial_conditions = get_initial_condition_for_event(
#'    puysegur_data, event_ID=250)
#' plot(initial_conditions)
#'
get_initial_condition_for_event<-function(source_zone_events_data, event_ID,
    force_file_download=FALSE){

    if(length(event_ID) != 1){
        stop('must have (length(event_ID)==1) in get_initial_condition_for_event -- you cannot pass a vector argument')
    }
    # Shorthand notation
    szed = source_zone_events_data

    # Get all raster names
    sz_rasters = szed$unit_source_statistics$initial_condition_file
   
    # Get event specific information 
    event_data = szed$events[event_ID,]
    event_slip = event_data$slip

    if(any(event_data$rate_annual == 0)){
        print('Warning: You requested the initial condition for an event that has an annual rate of zero (i.e. it is treated as impossible for the purposes of the PTHA)!')
    }

    event_raster_indices = scan(
        text=gsub("-", " ", event_data$event_index_string), quiet=TRUE)

    event_rasters = sz_rasters[event_raster_indices]

    # Figure out the raster names on this machine
    event_rasters_base = unlist(lapply(
        as.list(event_rasters), 
        f<-function(x) strsplit(x, split='SOURCE_ZONES')[[1]][2]))
    event_rasters_base = paste0('SOURCE_ZONES', event_rasters_base)

    if(!file.exists(event_rasters_base[1])){
        dir.create(dirname(event_rasters_base[1]), recursive=TRUE, showWarnings=FALSE)
    }

    # Figure out the raster names on NCI
    event_rasters_online = paste0(config_env$.GDATA_HTTP_BASE_LOCATION, 
        event_rasters_base)

    # Loop over all rasters, and if we can't find them in local directories,
    # then download them
    for(i in 1:length(event_rasters_base)){
        local_raster_name = event_rasters_base[i]
        if(force_file_download | (!file.exists(local_raster_name))){
            if(file.exists(event_rasters_online[i])){
                # This can happen if running code from NCI
                file.copy(event_rasters_online[i], local_raster_name)
            }else{
                download.file(event_rasters_online[i], local_raster_name)
            }
        }
    }

    # Read and sum the rasters
    variable_slip = ('event_slip_string' %in% names(event_data))
    if(variable_slip){
        slip = as.numeric(strsplit(event_data$event_slip_string, '_')[[1]])
    }else{
        slip = rep(event_data$slip, length(event_rasters_base))
    }

    r1 = raster(event_rasters_base[1])*0
    for(i in 1:length(event_rasters_base)){
        r1 = r1 + slip[i] * raster(event_rasters_base[i])
    }

    return(r1)
}

#' Download flow time-series from NCI
#'
#' @param source_zone_events_data output from \code{get_source_zone_events_data}
#' @param event_ID The row indices of the earthquake event(s) in source_zone_events_data$events
#' @param hazard_point_ID The numeric ID of the hazard point
#' @param target_polygon A SpatialPolygons object. All gauges inside this are selected
#' @param target_points A matrix of lon/lat point locations. The nearest gauge to each is selected
#' @param target_indices A vector with integer indices corresponding to where the gauge values are stored.
#' @param store_by_gauge Return the flow variables as a list with one gauge per
#' entry. Otherwise, return as a list with one event_ID per entry
#' @return Flow time-series
#'
get_flow_time_series_at_hazard_point<-function(source_zone_events_data, event_ID, 
    hazard_point_ID = NULL, target_polygon = NULL, target_points=NULL, target_indices = NULL,
    store_by_gauge=TRUE){

    is_null_hpID = is.null(hazard_point_ID)
    is_null_target_poly = is.null(target_polygon)
    is_null_target_points = is.null(target_points)
    is_null_target_indices = is.null(target_indices)

    only_one_input = is_null_hpID + is_null_target_poly + is_null_target_points + is_null_target_indices
    if(only_one_input != 3){
        print(only_one_input)
        stop('Only one of hazard_point_ID, target_polygon, target_points, target_indices should provided as non-NULL')
    }


    szed = source_zone_events_data

    if(!any(grepl('rptha', .packages(all=TRUE)))){
        stop('This function requires the rptha package to be installed, but the latter cannot be detected')
    }

    # Case of user-provided point IDs
    if(!is.null(hazard_point_ID)){
        # Find the index of the points matching event_ID inside the netcdf file
        indices_of_subset = get_netcdf_gauge_index_matching_ID(
            szed$gauge_netcdf_files[1],
            hazard_point_ID)

    }

    # Case of user-provided polygon
    if(!is.null(target_polygon)){

        indices_of_subset = get_netcdf_gauge_indices_in_polygon(
            szed$gauge_netcdf_files[1], target_polygon)

    }

    # Case of user-provided point locations
    if(!is.null(target_points)){
        indices_of_subset = get_netcdf_gauge_indices_near_points(
            szed$gauge_netcdf_files[1], target_points)
    }

    if(!is.null(target_indices)){
        indices_of_subset = target_indices
    }

    event_times = get_netcdf_gauge_output_times(szed$gauge_netcdf_files[1])
    gauge_locations = get_netcdf_gauge_locations(szed$gauge_netcdf_files[1], indices_of_subset)

    
    if(any(szed$events$rate_annual[event_ID] == 0)){
        print('Warning: You requested the initial condition for an event that has an annual rate of zero (i.e. it is treated as impossible for the purposes of the PTHA!)')
    }


    flow_var_batch = make_tsunami_event_from_unit_sources(
        szed$events[event_ID,], 
        szed$unit_source_statistics, 
        szed$gauge_netcdf_files,
        indices_of_subset=indices_of_subset, 
        verbose=FALSE,
        summary_function=NULL)

    if(store_by_gauge){
        # Store as a list, with one gauge in each list
        flow_var = vector(mode='list', length=length(indices_of_subset))
        for(i in 1:length(indices_of_subset)){
            flow_var[[i]] = array(0, 
                dim=c(length(event_ID), dim(flow_var_batch[[1]])[2:3]))
            for(j in 1:length(event_ID)){
                flow_var[[i]][j,,] = flow_var_batch[[j]][i,,]
            }
        }

        names(flow_var) = gauge_locations$gaugeID
    }else{
        # Store as a list, with one event in each entry
        names(flow_var_batch) = event_ID
        flow_var = flow_var_batch
    }

    output = list(time=event_times, flow=flow_var, locations=gauge_locations, events=szed$events[event_ID,])

    return(output)
}

#
# Utility function to determine the index of a hazard point in a netcdf file
#
# I often allow the user to pass a gaugeID, or a point, or an index. Then I need
# code to convert that into an index for the given netcdf file, and check that inputs are valid, etc. 
# This function takes care of that.
#
parse_ID_point_index_to_index<-function(netcdf_file, hazard_point_gaugeID, target_point, target_index){

    null_index = is.null(target_index)
    null_point = is.null(target_point)[1]
    null_ID = is.null(hazard_point_gaugeID)[1]

    only_one_input = null_ID + null_point + null_index
    if(only_one_input != 2){
        print(only_one_input)
        stop('Only one of hazard_point_gaugeID, target_point, target_index should provided as non-NULL')
    }

    if(null_index){
        if(!null_point){
            if(length(target_point) != 2){
                stop('Can only provide a single point')
            }
            # Find a site nearest the point
            target_index = get_netcdf_gauge_indices_near_points(
                netcdf_file,
                cbind(target_point[1], target_point[2]))
        }
    
        if(!null_ID){
            if(length(hazard_point_gaugeID) != 1){
                stop('Can only provide a single gauge ID')
            }
            # Find a site matching the ID
            target_index = get_netcdf_gauge_index_matching_ID(
                netcdf_file,
                hazard_point_gaugeID)
        }
    }

    return(target_index)
}


#'
#' Get the stage vs exceedance rate curves at a hazard point
#'
#' @param hazard_point_gaugeID numerical gaugeID of the hazard point of interest
#' @param target_point vector with c(lon, lat) of the target point
#' @param target_index integer index of the hazard point in the file
#' @param source_name Name of source-zone. If NULL, then return the rates for
#'     the sum over all source-zones
#' @param make_plot If TRUE, plot the stage vs return period curve for stochastic slip
#' @param non_stochastic_slip_sources If TRUE, also return curves for uniform
#'     and variable_uniform slip events
#' @param percentile_version either 'DG19' or 'DG18'. This affects the stage-percentile uncertainties. The 2019
#'     paper \url{https://doi.org/10.1007/s00024-019-02299-w} revised the method
#'     for computing percentile uncertainties in the maximum-stage. The latter results are
#'     used when version='DG19' (default). Otherwise one may use the results from the original
#'     PTHA18 report \url{http://dx.doi.org/10.11636/Record.2018.041} by setting version='DG18'.
#' @return list containing return period info for the source-zone
#'
get_stage_exceedance_rate_curve_at_hazard_point<-function(
    hazard_point_gaugeID = NULL, 
    target_point =NULL, 
    target_index = NULL,
    source_name = NULL,
    make_plot = FALSE,
    non_stochastic_slip_sources=FALSE,
    only_mean_rate_curve=FALSE,
    percentile_version = 'DG19'){

    # Depending on the percentile_version, specify a prefix which chooses the
    # correct file.
    if(percentile_version == 'DG19'){
        file_prefix = 'revised1_'
    }else if(percentile_version == 'DG18'){
        file_prefix = ''
    }else{
        stop('Version must be either "DG19" or "DG18"')
    }

    if(is.null(source_name)){
        source_name = 'sum_over_all_source_zones'
        stage_exceedance_rate_curves_file = paste0(config_env$.GDATA_OPENDAP_BASE_LOCATION,
            'EVENT_RATES/', file_prefix, 'tsunami_stage_exceedance_rates_', source_name, '.nc')
    }else{
        stage_exceedance_rate_curves_file = paste0(config_env$.GDATA_OPENDAP_BASE_LOCATION,
            'SOURCE_ZONES/', source_name, '/TSUNAMI_EVENTS/', file_prefix, 
            'tsunami_stage_exceedance_rates_', source_name, '.nc')
    }

    # Parse the input arguments into a target index
    target_index = parse_ID_point_index_to_index(
        stage_exceedance_rate_curves_file, hazard_point_gaugeID, 
        target_point, target_index)

    fid = nc_open(stage_exceedance_rate_curves_file, readunlim=FALSE)

    output = list()
    output$stage = fid$dim$stage$vals 

    vars = 'stochastic_slip_rate'
    if(non_stochastic_slip_sources){
        vars = c(vars, 'uniform_slip_rate', 'variable_uniform_slip_rate')
    }
    if(!only_mean_rate_curve){
        vars = c(vars, paste0(vars,'_upper_ci'), paste0(vars, '_lower_ci'), 
            paste0(vars, '_median'), paste0(vars, '_16pc'), paste0(vars, '_84pc'))
    }
    # Add in the variable_mu
    vars = c(vars, paste0('variable_mu_', vars))

    # Read the file
    for(i in 1:length(vars)){
        output[[vars[i]]] = ncvar_get(fid, vars[i],  start=c(1, target_index), count=c(-1,1))
    }

    if(!only_mean_rate_curve){
        ncdf_file_stations_only = get_file_with_gauges_only_if_on_NCI_THREDDS(stage_exceedance_rate_curves_file)
        fid2 = nc_open(ncdf_file_stations_only, readunlim=FALSE)
        output$lon = ncvar_get(fid2, 'lon', start=target_index, count=1)
        output$lat = ncvar_get(fid2, 'lat', start=target_index, count=1)
        output$elev = ncvar_get(fid2, 'elev', start=target_index, count=1)
        output$gaugeID = ncvar_get(fid2, 'gaugeID', start=target_index, count=1)
        nc_close(fid2)
    }
    output$target_index = target_index
    output$source_name = source_name
    output$stage_exceedance_rate_curves_file = stage_exceedance_rate_curves_file

    nc_close(fid)
    
    if(make_plot){
        title = paste0('Tsunami max-stage exceedance rates (stochastic slip, ', source_name, ')\n',
            'site = (', round(output$lon,3), ',', round(output$lat, 3), '); depth = ', 
            round(output$elev, 1), ' m; ID = ', round(output$gaugeID,3)) 
        options(scipen=5)
        plot(output$stage, output$stochastic_slip_rate, log='xy',
            xlim=c(0.02, 20), ylim=c(1.0e-04, 1), t='o', col='red', pch=19, cex=0.3,
            xlab='Maximum tsunami waterlevel above MSL', 
            ylab = 'Exceedance Rate (events/year)',
            main = title)
        points(output$stage, output$stochastic_slip_rate_upper, t='l', col='red')
        points(output$stage, output$stochastic_slip_rate_lower, t='l', col='red')
        points(output$stage, output$stochastic_slip_rate_median,t='l', col='red')
        points(output$stage, output$stochastic_slip_rate_16pc,  t='l', col='red')
        points(output$stage, output$stochastic_slip_rate_84pc,  t='l', col='red')

        points(output$stage, output$variable_mu_stochastic_slip_rate,       t='o', col='purple', cex=0.3, pch=19)
        points(output$stage, output$variable_mu_stochastic_slip_rate_upper, t='l', col='purple')
        points(output$stage, output$variable_mu_stochastic_slip_rate_lower, t='l', col='purple')
        points(output$stage, output$variable_mu_stochastic_slip_rate_median,t='l', col='purple')
        points(output$stage, output$variable_mu_stochastic_slip_rate_16pc,  t='l', col='purple')
        points(output$stage, output$variable_mu_stochastic_slip_rate_84pc,  t='l', col='purple')

        abline(v=c(0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20), col='orange', lty='dotted')
        abline(h=10**(seq(-6,0)), col='orange', lty='dotted')
        legend('topright', 
            c('Mean estimated rate', '2.5, 16, 50, 84, 97.5 % Percentiles', 
              'Mean estimated rate (variable mu)', '2.5, 16, 50, 84, 97.5 % Percentiles (variable mu)'
            ), 
            col=c('red', 'red', 'purple', 'purple'), 
            lty=c('solid', 'solid', 'solid', 'solid'), pch=c(19, NA, 19, NA), 
            pt.cex=c(0.3, NA, 0.3, NA), bg='white')
    }
    return(output)
}


#' Get stage vs exceedance rate curve for EVERY source-zone
#'
#' We loop over 'get_stage_exceedance_rate_curve_at_hazard_point'
#'
#' @param hazard_point_gaugeID numerical gaugeID of the hazard point of interest
#' @param target_point vector with c(lon, lat) of the target point
#' @param target_index integer index of the hazard point in the file
#' @param non_stochastic_slip_sources If TRUE, also get uniform and variable uniform rate curves
#' @param only_mean_rate_curve If TRUE, do not get credible interval rate curves
#' @return List of lists with the return period info for each source-zone
#'
get_stage_exceedance_rate_curves_all_sources<-function(
    hazard_point_gaugeID = NULL, 
    target_point = NULL, 
    target_index = NULL,
    non_stochastic_slip_sources=FALSE,
    only_mean_rate_curve=FALSE){

    all_sources = config_env$source_names_all

    outputs = vector(mode='list', length=length(all_sources))

    for(i in 1:length(all_sources)){

        if(i == 1){
            # On the first pass, we might not pass target_index
            outputs[[i]] = get_stage_exceedance_rate_curve_at_hazard_point(
                hazard_point_gaugeID,
                target_point,
                target_index, 
                source_name = all_sources[i],
                non_stochastic_slip_sources = non_stochastic_slip_sources,
                only_mean_rate_curve=only_mean_rate_curve)
        }else{
            # On the second pass, we can pass target_index, which is faster
            outputs[[i]] = get_stage_exceedance_rate_curve_at_hazard_point(
                target_index = outputs[[1]]$target_index, 
                source_name = all_sources[i],
                non_stochastic_slip_sources = non_stochastic_slip_sources,
                only_mean_rate_curve=only_mean_rate_curve)
        }
        names(outputs)[i] = all_sources[i]

    }
    return(outputs)
}

#'
#' For a point, get a list which contains (for each source-zone),
#' peak_stage, and (optionally) Mw, variable_mu_Mw, and information on
#' whether the scenario rate is nonzero. This can be helpful to 
#' select events for further study.
#'
#' @param hazard_point_gaugeID numerical gaugeID of the hazard point of interest
#' @param target_point vector with c(lon, lat) of the target point
#' @param target_index integer index of the hazard point in the file
#' @param all_source_names vector with the source names to extract data from
#' @param slip_type 'stochastic' or 'variable_uniform' or 'uniform'
#' @param include_earthquake_data TRUE/FALSE do we also read some earthquake information?
#' @param max_tries if a download fails, try again this many times. This can help when
#' using unreliable internet connections and/or servers.
#' @return list (one entry for each source) containing a list with 
#'   peak_stage, and (if include_earthquake_data=TRUE) Mw, event_rate, ...
#'
get_peak_stage_at_point_for_each_event<-function(hazard_point_gaugeID = NULL, 
    target_point=NULL, target_index = NULL, all_source_names = NULL,
    slip_type = 'stochastic', include_earthquake_data=TRUE, max_tries=20){


    # Here it doesn't matter whether we use the revised1_XXXX file, or the original file
    # that doesn't have the 'revised1_' prefix, because the data we draw is the same.
    stage_exceedance_rate_curves_file = paste0(config_env$.GDATA_OPENDAP_BASE_LOCATION,
        'EVENT_RATES/revised1_tsunami_stage_exceedance_rates_sum_over_all_source_zones.nc')
    # Parse the input arguments into a target index
    target_index = parse_ID_point_index_to_index(
        stage_exceedance_rate_curves_file, hazard_point_gaugeID, 
        target_point, target_index)

    if(is.null(all_source_names)){
        #all_source_names = basename(dirname(Sys.glob('SOURCE_ZONES/*/EQ_SOURCE')))
        all_source_names = config_env$source_names_all
    }

    if(slip_type == 'stochastic'){
        file_base = 'all_stochastic_slip_earthquake_events_'
    }else if(slip_type == 'variable_uniform'){
        file_base = 'all_variable_uniform_slip_earthquake_events_'
    }else if(slip_type == 'uniform'){
        file_base = 'all_uniform_slip_earthquake_events_'
    }else{
        stop('unrecognized slip_type')
    }

    output = vector(mode='list', length=length(all_source_names))
    names(output) = all_source_names

    # Download the data for all sourcezones
    for(i in 1:length(all_source_names)){
    
        nm = all_source_names[i]
        print(nm)
        try_again = TRUE

        # Allow for the download to fail up to 'max_tries' times
        counter = 0
        #max_tries = 20
        has_vars = c(FALSE, !include_earthquake_data)
        while(try_again){

            counter = counter + 1

            # Read the max stage
            if(!has_vars[1]){
                # Read a new version of the files containing ony 'max_stage', for speed.
                nc_file1 = paste0(config_env$.GDATA_OPENDAP_BASE_LOCATION, 'SOURCE_ZONES/',
                    nm, '/TSUNAMI_EVENTS/MAX_STAGE_ONLY_', file_base, 'tsunami_', 
                    nm, '_MAX_STAGE_ONLY.nc') 
                fid1 = nc_open(nc_file1, readunlim=FALSE)
                local_max_stage = try(ncvar_get(fid1, 'max_stage', start=c(1,target_index), 
                    count=c(fid1$dim$event$len,1)))
                nc_close(fid1)
                # Record if it failed, to try again later
                #if((class(local_max_stage) != 'try-error')) has_vars[1] = TRUE
                if(!is(local_max_stage, 'try-error')) has_vars[1] = TRUE
            }

            # Read Mw and the event rate from the file that doesn't contain the tsunami
            # max-stage, because the read access is faster
            if(!has_vars[2]){
                nc_file2 = paste0(config_env$.GDATA_OPENDAP_BASE_LOCATION, 'SOURCE_ZONES/',
                    nm, '/TSUNAMI_EVENTS/', file_base, 
                    nm, '.nc')
                fid2 = nc_open(nc_file2, readunlim=FALSE, suppress_dimvals=TRUE)
                local_Mw = try(ncvar_get(fid2, 'Mw'))
                local_rate = try(ncvar_get(fid2, 'rate_annual'))
                local_Mw_variable_mu = try(ncvar_get(fid2, 'variable_mu_Mw'))
                local_rate_variable_mu = try(ncvar_get(fid2, 'variable_mu_rate_annual'))
                nc_close(fid2)
                # Record if it failed, to try again later
                #if((class(local_Mw) != 'try-error') & 
                #   (class(local_rate) != 'try-error') &
                #   (class(local_rate_variable_mu) != 'try-error') &
                #   (class(local_Mw_variable_mu) != 'try-error')
                #    ) has_vars[2] = TRUE
                if((!is(local_Mw, 'try-error')) & 
                   (!is(local_rate, 'try-error')) &
                   (!is(local_rate_variable_mu, 'try-error')) &
                   (!is(local_Mw_variable_mu, 'try-error'))
                    ) has_vars[2] = TRUE

            }

            if(include_earthquake_data){
                #
                # Full output case
                #
                output[[i]] = list(
                    Mw = local_Mw,
                    max_stage = local_max_stage,
                    #period = local_period,
                    scenario_rate_is_positive = (local_rate>0),
                    target_index=target_index,
                    slip_type=slip_type,
                    variable_mu_Mw = local_Mw_variable_mu,
                    variable_mu_scenario_rate_is_positive = (local_rate_variable_mu>0)
                    )

            }else{
                #
                # Minimal output case
                #

                output[[i]] = list(
                    max_stage = local_max_stage,
                    target_index=target_index,
                    slip_type=slip_type
                    )

            }


            # Error handling
            if(!all(has_vars == TRUE)){

                if(counter <= max_tries){
                    try_again = TRUE
                    print('remote read failed, trying again')
                }

                if(counter > max_tries){
                    try_again = FALSE
                    print('remote read failed too many times, skipping')
                }
            }else{

                try_again = FALSE
            }
        } # End of this source zone
    } # End loop over all source-zones

    return(output)
}

#'
#' Function to summarize event properties. This MAY help
#' in choosing one or a few scenarios from a set of scenarios. 
#'
summarise_events<-function(events_near_desired_stage){

    # shorthand
    ends = events_near_desired_stage
   
    mws = ends$events$Mw
    rates = ends$events$rate_annual
    peak_slip = sapply(ends$events$event_slip_string, 
        f<-function(x) max(as.numeric(strsplit(x, '_')[[1]])),
        USE.NAMES=FALSE)
    mean_slip = sapply(ends$events$event_slip_string, 
        f<-function(x) mean(as.numeric(strsplit(x, '_')[[1]])),
        USE.NAMES=FALSE)
    nsources = sapply(ends$events$event_slip_string, 
        f<-function(x) length(as.numeric(strsplit(x, '_')[[1]])),
        USE.NAMES=FALSE)
    peak_slip_alongstrike = ends$events$peak_slip_alongstrike_ind 

    magnitude_prop_le = sapply(mws, f<-function(x) sum(rates * (mws <= x)))/sum(rates)
    magnitude_prop_lt = sapply(mws, f<-function(x) sum(rates * (mws < x)))/sum(rates)

    magnitude_prop_mid = 0.5*(magnitude_prop_le + magnitude_prop_lt)

    # Reproducible jitter
    if(exists('.Random.seed')){
        oldseed = .Random.seed
    }else{
        p = runif(1) # Now we will have a random seed
        oldseed = .Random.seed
    }
    set.seed(1)
 
    obs = cbind(jitter(magnitude_prop_mid, 0.0001), peak_slip, mean_slip, nsources, peak_slip_alongstrike)
    if(nrow(obs) > 1){
        # We can do a mahalanobis distance calculation
        obs = apply(obs, 2, f<-function(x) qnorm(rank(x)/(length(x)+1)))
        mh_distance = try(mahalanobis(obs, center=mean(obs), cov=cov(obs)), silent=TRUE)
        #if(class(mh_distance) == 'try-error'){
        if(is(mh_distance, 'try-error')){
            # Mahalanobis distance failed (perhaps too few events, or
            # problematic input data leading to singular covariance matrix).
            # Prioritize events with intermediate magnitudes
            mh_distance = abs(magnitude_prop_mid - 0.5) 
        }
    }else{
        mh_distance = 0
    }

    # Undo reproducible random jitter
    .Random.seed <<- oldseed

    return(data.frame(mws, peak_slip, mean_slip, nsources, peak_slip_alongstrike, 
        magnitude_prop_mid, mh_distance))
}

#'
#' Download a DEM with 1 in 'below MSL' and 0 in 'above MSL' regions. This
#' was derived from the input merged DEM used for the PTHA18, created with
#' this script (which in turn relies on the other codes in the same directory): 
#'
#' https://github.com/GeoscienceAustralia/ptha/blob/master/R/examples/austptha_template/DATA/ELEV/merged_dem/make_wet_or_dry_dem.R
#' 
get_wet_or_dry_DEM<-function(force_download_again=FALSE){
    
    wet_or_dry_DEM_file = paste0(config_env$.GDATA_HTTP_BASE_LOCATION, 
        'DATA/wet_or_dry_gebco_ga250_dem_patched.tif')

    output_file = './.wet_or_dry_gebco_ga250_dem_patched/wet_or_dry_gebco_ga250_dem_patched.tif'

    if(file.exists(output_file) & !force_download_again){
        # We do not need to download the data
        wd = raster(output_file)
    }else{
        # We do need to download the data
        dir.create(dirname(output_file), showWarnings=FALSE)
        if(file.exists(wet_or_dry_DEM_file)){
            file.copy(wet_or_dry_DEM_file, output_file, overwrite=force_download_again)
        }else{
            download.file(wet_or_dry_DEM_file, output_file)
        }
        wd = raster(output_file)
        writeRaster(wd, file=output_file, options=c('COMPRESS=DEFLATE'), 
                    overwrite=force_download_again)
    }


    return(wd)
}

#' Read the initial potential energy from the precomputed database
#'
#' For details on the potential energy calculation, see here: https://github.com/GeoscienceAustralia/ptha/blob/master/R/examples/austptha_template/SOURCE_ZONES/compute_initial_condition_potential_energy.R
#'
#' @param source_zone Name of source_zone
#' @param slip_type 'stochastic' or 'variable_uniform' or 'uniform'
#' @param desired_event_rows integer vector giving the rows of the table that
#' are desired. If NULL, read all rows (unless range_list is not NULL, see below)
#' @param chunk_size The chunk_size passed to read_table_from_netcdf. This can impact the efficiency when only
#' reading a subset of the file. A small chunk size will lead to many reads, whereas a large chunk size may lead
#' to too much data being read at once. The best size depends on the system and the interconnect.
#'
get_source_zone_events_potential_energy<-function(source_zone, slip_type='stochastic', 
    desired_event_rows=NULL, chunk_size=1000){

    # First check that a valid source-zone was provided
    err = FALSE
    if(is.null(source_zone)){
        err = TRUE
    }else{
        if(sum(config_env$source_names_all == source_zone) == 0) err=TRUE
    }

    if(err){
        print('You did not pass a valid source_zone to get_source_zone_events_data. The allowed source_zone values are:')
        print(paste0('   ', config_env$source_names_all))
        print('Please pass one of the above source_zone names to this function to get its metadata')
        # Fail gracefully
        output = list(initial_potential_energy = NA, rate_annual=NA, 
                      weight_with_nonzero_rate = NA, Mw = NA)
        return(invisible(output))
    }

    library(rptha)

    nc_web_addr = paste0(config_env$.GDATA_OPENDAP_BASE_LOCATION, 
        'SOURCE_ZONES/', source_zone, '/TSUNAMI_EVENTS/POTENTIAL_ENERGY_all_', slip_type, 
        '_slip_earthquake_events_', source_zone, '_POTENTIAL_ENERGY.nc')

    potential_energy_data = read_table_from_netcdf(nc_web_addr, 
        desired_rows=desired_event_rows, chunk_size=chunk_size)

    return(potential_energy_data)

}

#' Randomly sample scenarios given their magnitudes and rates
#'
#' Generate a random sample of PTHA18 scenarios from a source-zone, with
#' sampling stratified by magnitude. Within each magnitude, the chance of each
#' event being sampled is proportional to the event_rate, optionally multiplied
#' by a user-specified event_importance (defaults to 1). 
#'
#' @param event_rates vector of PTHA18 scenario rates (one for each scenario)
#' @param event_Mw vector of PTHA18 scenario magnitudes (one for each scenario) 
#' The values must be in 7.2, 7.3, ...9.6, 9.7, 9.8 as for PTHA18.
#' @param event_importance if not NULL, then a vector of non-negative numbers
#' (one for each scenario) giving the "importance" of each scenario. By default
#' this is treated as a constant (=1). For each magnitude, the probability of
#' sampling the scenario is proportional to "event_rates * event_importance",
#' so events with higher importance are more likely to be sampled than they
#' otherwise would be. The output scenario rates are adjusted to account for
#' this, so the resulting set of scenarios and rates can still be treated as a
#' random sample from the source-zone. But this sample will better resolve
#' scenarios with higher importance, and more poorly resolve scenarios with
#' lower importance.
#' @param samples_per_Mw a function returning the number of scenarios to sample
#' for each magnitude.
#' @param mw_limits only sample from magnitudes greater than mw_limits[1], and
#' less than mw_limits[2]
#' @param return_as_table convert the output to a data.frame, rather than a list of lists.
#' @return a data.frame (if return_as_table=TRUE) or a list of lists (one per magnitude). 
#' Each contains the overall rate of scenarios with the given magnitude
#' (rate_with_this_mw), the magnitude (mw), the random scenario indices
#' corresponding to event_rates and event_Mw (inds), the product of the
#' rate_with_this_mw and the standard importance sampling weights for each
#' scenario (importance_sampling_scenario_rates), the product of the
#' rate_with_this_mw and the self-normalised importance sampling weights for
#' each scenario (importance_sampling_scenario_rates_self_normalised), the
#' standard importance sampling weights for each scenario, which may not sum to
#' 1, but do so on average, and lead to unbiased integral estimates
#' (importance_sampling_scenario_weights), and the self-normalised importance
#' sampling weights for each scenario, which will sum to 1, but lead to
#' integral estimates that can be biased (but are asymptotically unbiased).
#' @details The function returns either a data.frame or a list of lists (one per unique magnitude)
#' containing the random scenario indices, and associated rates that can be
#' assigned to each scenario for consistency with the PTHA18. The rates are
#' computed with importance sampling, and provide statistical consistency with
#' the PTHA18 even if a non-uniform event_importance is specified. Two
#' alternative importance-sampling based rates are provided. These only differ
#' when the event_importance is specified and non-uniform within a magnitude.
#' One partitions the rate using "standard importance sampling weights":
#'     (1/number_of_random_scenarios) * [ 
#'       ( PTHA18 conditional probability for the magnitude ) / 
#'       ( (event_rates * event_importance) for the magnitude, normalised to a density)
#'         ]  
#' where the terms in [ ] are evaluated at the randomly selected scenarios. The
#' other uses "self-normalised importance sampling weights":
#'     [ (1/event_importance_for_the_randomly_selected_scenarios )/
#'     sum( 1/event_importance_for_the_randomly_selected_scenarios) ]
#' To partition the rates among scenarios for each Mw, we multiply these by the
#' rate of scenarios with the given Mw. Each method has different strengths and
#' weaknesses. The former importance sampling based rates can be used to
#' compute unbiased estimates of integrals, whereas the latter is only
#' asymptotically unbiased as the number of scenarios approaches infinity.
#' However the standard importance sampling weights do not necessarily sum to
#' 1, so lead to inconsistencies with the PTHA18 rates for each magnitude,
#' whereas the self-normalised importance sampling weights retain consistency
#' with the PTHA18 rates for each magnitude. Depending on the application, one
#' may prefer one or the other.
#' 
randomly_sample_scenarios_by_Mw_and_rate<-function(
    event_rates,
    event_Mw,
    event_importance = NULL,
    samples_per_Mw=function(Mw){round(50 + 0*Mw)},
    mw_limits=c(7.15, 9.85),
    return_as_table=TRUE){

    unique_Mw = sort(unique(event_Mw))

    diff_unique_Mw = diff(unique_Mw)

    if(length(diff_unique_Mw) > 1){

        # Confirm that our expectations of PTHA18 scenarios are met.
        if(any(abs(diff_unique_Mw - diff_unique_Mw[1]) > 1.0e-06) ){
            stop('event_Mw is not evenly spaced')
        }
    }

    # Ignore small scenarios, and scenarios exceeding Mw-max
    unique_Mw = unique_Mw[ ((unique_Mw > mw_limits[1]) & 
                            (unique_Mw < mw_limits[2])) ] 

    # By default give all scenarios equal importance.
    if(is.null(event_importance)) event_importance = event_Mw * 0 + 1

    # Sanity check on inputs
    if( any( (event_rates < 0) | (event_importance < 0)) ){
        stop('event_rates and event_importance must be nonnegative')
    }

    random_scenario_info = lapply(unique_Mw,
        f<-function(mw){
            # Match Mw [careful with real numbers]
            k = which(abs(event_Mw - mw) < 1.0e-03)
            if(length(k) <= 1){
                stop('error: Only one scenario -- need to be careful using sample below')
            }

            # Get the rate of any event with this mw
            rate_with_this_mw = sum(event_rates[k])
            nsam = round(samples_per_Mw(mw))

            if((nsam > 0) & sum(event_rates[k] * event_importance[k]) > 0){

                local_sample = sample(1:length(k), size=nsam, 
                    prob=event_rates[k]*event_importance[k], 
                    replace=TRUE)
                sample_of_k = k[local_sample]

                # The original scenario conditional probability distribution
                dist_f = event_rates[k]/sum(event_rates[k])
                # The distribution we sampled from above
                dist_g = (event_rates[k]*event_importance[k]) / 
                      sum(event_rates[k]*event_importance[k])
                # The exact importance-sampling correction -- while these
                  # weights do not sum to 1, they make the estimators unbiased
                alternate_weights = (1/length(local_sample)) *
                    (dist_f[local_sample] / dist_g[local_sample])

                # Importance sampling correction -- this is the
                # 'self-normalised' approach.
                # Estimate a rate for each individual scenario such that the
                # sampled scenario weights add to the target weight, and
                # reflect the original distribution
                # Their chance of being sampled was inflated by
                # 'event_importance', so we deflate by 'event_importance' here.
                random_scenario_rates = rate_with_this_mw * 
                    (1/event_importance[sample_of_k]) / 
                    sum((1/event_importance[sample_of_k]))

                # Importance sampling correction -- although "alternate_weights"
                # does not sum to 1, the self normalised importance sampling
                # leads to some finite-sample-bias in integral estimates,
                # whereas this approach is unbiased even in finite samples.
                alternate_random_scenario_rates = 
                    rate_with_this_mw * alternate_weights

                return(data.frame(
                    inds = sample_of_k, 
                    mw = rep(mw, nsam),
                    rate_with_this_mw = rep(rate_with_this_mw, nsam),
                    # Regular importance sampling
                    importance_sampling_scenario_rates = 
                        alternate_random_scenario_rates, 
                    # Self-normalised importance sampling
                    importance_sampling_scenario_rates_self_normalised = 
                        random_scenario_rates, 
                    importance_sampling_scenario_weights = alternate_weights,
                    importance_sampling_scenario_weights_self_normalised = 
                        random_scenario_rates/rate_with_this_mw 
                    ))
            }else{
                # Case with all rates=0
                return(data.frame(
                    inds = NA, 
                    mw = mw,
                    rate_with_this_mw = rate_with_this_mw,
                    importance_sampling_scenario_rates = NA, 
                    importance_sampling_scenario_rates_self_normalised = NA,
                    importance_sampling_scenario_weights = NA,
                    importance_sampling_scenario_weights_self_normalised = NA
                    ))
            }
        })

    if(return_as_table){
        if(length(random_scenario_info) > 1){
            random_scenario_info = do.call(rbind, random_scenario_info)
        }
    }

    return(random_scenario_info)
}

