suppressPackageStartupMessages(library(raster))
if(!exists('.HAVE_SOURCED_CONFIG')) source('R/config.R', local=TRUE, chdir=FALSE)

#' Read key summary statistics for earthquake events on the source-zone
#'
#' The code reads the files from the web, so an internet connection is required.
#' At the moment no caching of data is implemented.
#'
#' @param source_zone Name of source_zone
#' @return list with 'events' giving summary statistics for the earthquake
#' events, and 'unit_source_statistics' giving summary statistics for each
#' unit source
#' @export
#' @examples
#' puysegur_data = get_source_zone_events_data('puysegur')
#'
get_source_zone_events_data<-function(source_zone){

    #
    # Get the earthquake events data
    #
    csv_web_addr = paste0(.GDATA_HTTP_BASE_LOCATION, 'SOURCE_ZONES/', 
        source_zone, '/TSUNAMI_EVENTS/all_eq_events_', source_zone, '.csv')    

    events_data = read.csv(csv_web_addr, stringsAsFactors=FALSE)

    #
    # Get the unit source summary statistics
    #
    csv_web_addr = paste0(.GDATA_HTTP_BASE_LOCATION, 'SOURCE_ZONES/', 
        source_zone, '/TSUNAMI_EVENTS/unit_source_statistics_', source_zone, 
        '.csv')    

    unit_source_statistics = read.csv(csv_web_addr, stringsAsFactors=FALSE)

    #
    # Get the netcdf gauge files, by reading a text file that contains their
    # locations. We then have to re-order the result to match the row-order
    # in unit_source_statistics
    #
    text_web_addr = paste0(.GDATA_HTTP_BASE_LOCATION, 'SOURCE_ZONES/', 
        source_zone, '/TSUNAMI_UNIT_SOURCE/gauge_files_list.txt')    

    gauge_netcdf_files = readLines(text_web_addr)

    # Order the gauge files in the same way as the unit sources
    unit_source_flag = gsub('.tif', '/', basename(unit_source_statistics$initial_condition_file))
    gauge_netcdf_files_reordered = rep(NA, length(gauge_netcdf_files))
    for(i in 1:length(unit_source_flag)){
        # Find the netcdf file corresponding to the unit source
        ind = grep(unit_source_flag[i], gauge_netcdf_files)
        if(length(ind)!= 1) stop('Error in matching of gauge files and unit sources')

        gauge_netcdf_files_reordered[i] = gauge_netcdf_files[ind]
    }
  
    # Append the web address to the files 
    gauge_netcdf_files_reordered = paste0(.GDATA_OPENDAP_BASE_LOCATION, 'SOURCE_ZONES/', 
        source_zone, '/TSUNAMI_UNIT_SOURCE/', gauge_netcdf_files_reordered) 

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
    r1 = raster(event_rasters_base[1])
    if(length(event_rasters_base) > 1){
        for(i in 2:length(event_rasters_base)){
            r1 = r1 + raster(event_rasters_base[i])
        }
    }
    r1 = r1 * event_slip

    return(r1)
}

#' Download flow time-series from NCI
#'
#' @param source_zone_events_data output from \code{get_source_zone_events_data}
#' @param event_ID The row index of the earthquake event in source_zone_events_data$events
#' @param hazard_point_ID The numeric ID of the hazard point
#' @param target_polygon A SpatialPolygons object. All gauges inside this are selected
#' @param target_points A matrix of lon/lat point locations. The nearest gauge to each is selected
#' @param unpack_to_list Return flow_var as a list with one gauge per element
#' @return Flow time-series
get_flow_time_series_at_hazard_point<-function(source_zone_events_data, event_ID, 
    hazard_point_ID = NULL, target_polygon = NULL, target_points=NULL,
    unpack_to_list=TRUE){

    szed = source_zone_events_data

    if(!any(grepl('rptha', .packages(all=TRUE)))){
        stop('This function requires the rptha package to be installed, but the latter cannot be detected')
    }else{
        source('R/sum_tsunami_unit_sources.R', local=TRUE)
    }

    # Case of user-provided point IDs
    if(!is.null(hazard_point_ID)){
        if(!is.null(target_polygon) | !is.null(target_points)){
            stop('Only one of hazard_point_ID, target_polygon, target_points should provided as non-NULL')
        }

        # Find the index of the points matching event_ID inside the netcdf file
        event_indices = get_netcdf_gauge_index_matching_ID(
            szed$gauge_netcdf_files[1],
            hazard_point_ID)

    }

    # Case of user-provided polygon
    if(!is.null(target_polygon)){

        if(!is.null(hazard_point_ID) | !is.null(target_points)){
            stop('Only one of hazard_point_ID, target_polygon, target_points should provided as non-NULL')
        }
        
        event_indices = get_netcdf_gauge_indices_in_polygon(
            szed$gauge_netcdf_files[1], target_polygon)

    }

    # Case of user-provided point locations
    if(!is.null(target_points)){
        if(!is.null(hazard_point_ID) | !is.null(target_polygon)){
            stop('Only one of hazard_point_ID, target_polygon, target_points should provided as non-NULL')
        }
        event_indices = get_netcdf_gauge_indices_near_points(
            szed$gauge_netcdf_files[1], target_points)

    }

    event_times = get_netcdf_gauge_output_times(szed$gauge_netcdf_files[1])
    gauge_locations = get_netcdf_gauge_locations(szed$gauge_netcdf_files[1], event_indices)

    flow_var_batch = make_tsunami_event_from_unit_sources(
        szed$events[event_ID,], 
        szed$unit_source_statistics, 
        szed$gauge_netcdf_files,
        #get_flow_time_series_function = get_flow_time_series_SWALS,  
        indices_of_subset=event_indices, 
        verbose=FALSE,
        summary_function=NULL)[[1]]

    # Optionally output flow_var as a list
    if(unpack_to_list){
        flow_var = list()
        for(i in 1:dim(flow_var_batch)[1]){
            flow_var[[i]] = flow_var_batch[i,,,drop=FALSE]
        }
        names(flow_var) = gauge_locations$gaugeID
    }else{
        # Save the copy and keep in compact matrix format
        flow_var = flow_var_batch
    }

    output = list(time=event_times, flow=flow_var, locations=gauge_locations)

    return(output)
}
