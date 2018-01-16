#
# Codes to remotely access various datasets from new Australian PTHA, on the NCI
#

suppressPackageStartupMessages(library(raster))
if(!exists('.HAVE_SOURCED_CONFIG')) source('R/config.R', local=TRUE, chdir=FALSE)
source('R/sum_tsunami_unit_sources.R')

#' Read key summary statistics for earthquake events on the source-zone
#'
#' The code reads the files from the web, so an internet connection is required.
#' At the moment no caching of data is implemented. Also, the file sizes range
#' greatly (from a few MB up to 15GB). Subsetting must be used on the large files
#' (assuming typical 2018 internet speeds).
#'
#' @param source_zone Name of source_zone
#' @return list with 'events' giving summary statistics for the earthquake
#' events, and 'unit_source_statistics' giving summary statistics for each
#' unit source
#' @export
#' @examples
#' puysegur_data = get_source_zone_events_data('puysegur')
#'
get_source_zone_events_data<-function(source_zone, slip_type='uniform'){

    stopifnot(slip_type %in% c('uniform', 'stochastic', 'variable_uniform'))

    #
    # Get the earthquake events data
    #
    nc_web_addr = paste0(.GDATA_OPENDAP_BASE_LOCATION, 'SOURCE_ZONES/', 
        source_zone, '/TSUNAMI_EVENTS/all_', slip_type, 
        '_slip_earthquake_events_', source_zone, '.nc')    

    events_data = read_table_from_netcdf(nc_web_addr)

    #
    # Get the unit source summary statistics
    #
    nc_web_addr = paste0(.GDATA_OPENDAP_BASE_LOCATION, 'SOURCE_ZONES/', 
        source_zone, '/TSUNAMI_EVENTS/unit_source_statistics_', source_zone, 
        '.nc')    

    unit_source_statistics = read_table_from_netcdf(nc_web_addr)

    gauge_netcdf_files = unit_source_statistics$tide_gauge_file 

    # Append the web address to the files 
    gauge_netcdf_files_reordered = sapply(gauge_netcdf_files, adjust_path_to_gdata_base_location)

    output = list()
    output[['events']] = events_data
    output[['unit_source_statistics']] = unit_source_statistics
    output[['gauge_netcdf_files']] = gauge_netcdf_files_reordered
    

    return(invisible(output))
}

#' Create initial conditions for tsunami model
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

    # Shorthand notation
    szed = source_zone_events_data

    # Get all raster names
    sz_rasters = szed$unit_source_statistics$initial_condition_file
   
    # Get event specific information 
    event_data = szed$events[event_ID,]
    event_slip = event_data$slip

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
    event_rasters_online = paste0(.GDATA_HTTP_BASE_LOCATION, 
        event_rasters_base)

    # Loop over all rasters, and if we can't find them in local directories,
    # then download them
    for(i in 1:length(event_rasters_base)){
        local_raster_name = event_rasters_base[i]
        if(force_file_download | (!file.exists(local_raster_name))){
            download.file(event_rasters_online[i], local_raster_name)
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
    }else{
        source('R/sum_tsunami_unit_sources.R', local=TRUE)
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

    flow_var_batch = make_tsunami_event_from_unit_sources(
        szed$events[event_ID,], 
        szed$unit_source_statistics, 
        szed$gauge_netcdf_files,
        #get_flow_time_series_function = get_flow_time_series_SWALS,  
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
#' @param
get_stage_exceedance_rate_curve_at_hazard_point<-function(hazard_point_gaugeID = NULL, 
    target_point=NULL, target_index = NULL){

    stage_exceedance_rate_curves_file = paste0(.GDATA_OPENDAP_BASE_LOCATION,
        'EVENT_RATES/tsunami_stage_exceedance_rates_sum_over_all_source_zones.nc')

    # Parse the input arguments into a target index
    target_index = parse_ID_point_index_to_index(
        stage_exceedance_rate_curves_file, hazard_point_gaugeID, 
        target_point, target_index)

    fid = nc_open(stage_exceedance_rate_curves_file, readunlim=FALSE)

    output = list()
    output$stage = fid$dim$stage$vals 
    vars = c('stochastic_slip_rate', 'stochastic_slip_rate_upper_ci', 'stochastic_slip_rate_lower_ci',
        'uniform_slip_rate', 'uniform_slip_rate_upper_ci', 'uniform_slip_rate_lower_ci',
        'variable_uniform_slip_rate', 'variable_uniform_slip_rate_upper_ci', 'variable_uniform_slip_rate_lower_ci')
    # Read the file
    for(i in 1:length(vars)){
        output[[vars[i]]] = ncvar_get(fid, vars[i],  start=c(1, target_index), count=c(-1,1))
    }

    output$lon = ncvar_get(fid, 'lon', start=target_index, count=1)
    output$lat = ncvar_get(fid, 'lat', start=target_index, count=1)
    output$elev = ncvar_get(fid, 'elev', start=target_index, count=1)
    output$gaugeID = ncvar_get(fid, 'gaugeID', start=target_index, count=1)

    nc_close(fid)

    return(output)
}

#
# For a point, get a list which contains (for each source-zone),
# Mw, max_stage
#
get_peak_stage_at_point_for_each_event<-function(hazard_point_gaugeID = NULL, 
    target_point=NULL, target_index = NULL, all_source_names = NULL){


    stage_exceedance_rate_curves_file = paste0(.GDATA_OPENDAP_BASE_LOCATION,
        'EVENT_RATES/tsunami_stage_exceedance_rates_sum_over_all_source_zones.nc')
    # Parse the input arguments into a target index
    target_index = parse_ID_point_index_to_index(
        stage_exceedance_rate_curves_file, hazard_point_gaugeID, 
        target_point, target_index)

    if(is.null(all_source_names)){
        all_source_names = basename(dirname(Sys.glob('SOURCE_ZONES/*/EQ_SOURCE')))
    }

    output = vector(mode='list', length=length(all_source_names))
    names(output) == all_source_names

    for(i in 1:length(all_source_names)){
    
        nm = all_source_names[i]
        print(nm)

        nc_file1 = paste0(.GDATA_OPENDAP_BASE_LOCATION, 'SOURCE_ZONES/',
            nm, '/TSUNAMI_EVENTS/all_stochastic_slip_earthquake_events_tsunami_', 
            nm, '.nc')
        fid1 = nc_open(nc_file1, readunlim=FALSE)
        local_max_stage = ncvar_get(fid1, 'max_stage', start=c(1,target_index), 
            count=c(fid1$dim$event$len,1))
        nc_close(fid1)

        nc_file2 = paste0(.GDATA_OPENDAP_BASE_LOCATION, 'SOURCE_ZONES/',
            nm, '/TSUNAMI_EVENTS/all_stochastic_slip_earthquake_events_', 
            nm, '.nc')
        fid2 = nc_open(nc_file2, readunlim=FALSE)
        local_Mw = ncvar_get(fid2, 'Mw')
        nc_close(fid2)

        output[[i]] = list(
            Mw = local_Mw,
            max_stage = local_max_stage,
            target_index=target_index
            )

    }

    return(output)
}
