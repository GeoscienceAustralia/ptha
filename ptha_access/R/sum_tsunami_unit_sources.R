#
# Code for summing unit-source tsunami waveforms
#


suppressPackageStartupMessages(library(ncdf4))
suppressPackageStartupMessages(library(rptha))

# Get the config data -- this will allow us to use information on the NCI THREDDS filesystem
config_env = new.env()
source('R/config.R', local=config_env)

#'
#' Get the 'input_stage_raster' attribute from the tide gauge netcdf files. This
#' helps us match netcdf files to unit sources
#'
#' @param netcdf_file name of netcdf tide gauge file
#' @return value of the input_stage_raster attribute of the netcdf file
#' @export
get_netcdf_attribute_initial_stage_raster<-function(netcdf_file){

    fid = nc_open(netcdf_file, readunlim=FALSE)
    myatt = ncatt_get(fid, varid=0, 'input_stage_raster')$value
    nc_close(fid)
    return(myatt)
}

#' Undo precision issues in gaugeID (caused by storing as float)
#'
#' @param gaugeID vector of numeric gauge IDs, which will be close
#' to a number having 2 non-zero decimal places, but may be slightly
#' rounded due to floating point precision issues
#'
clean_gaugeID<-function(gaugeID){

    clean_gauge_val = round(gaugeID, 2)
    if(any(abs(clean_gauge_val - gaugeID) > 0.005)){
        stop('Error with gauge IDs')
    }
    return(clean_gauge_val)
}

#' It can be slow to read gauge locations remotely if using files where station
#' is an unlimited dimension. To work around this, the current function takes
#' a netcdf file, and if it deems the file is on NCI THREDDS then it returns the
#' name of another file which only contains the station lon/lat/elev/gaugeID.
#' Calling this in the right place can speed up the program
get_file_with_gauges_only_if_on_NCI_THREDDS<-function(netcdf_file){

    on_NCI = (length(grep(config_env$.GDATA_OPENDAP_BASE_LOCATION, netcdf_file, fixed=TRUE)) > 0)
    if(on_NCI){
        # Read a file that ONLY contains gauges -- this is generally faster
        netcdf_file = config_env$adjust_path_to_gdata_base_location(
            paste0(config_env$.GDATA_OPENDAP_BASE_LOCATION, 
                   'EVENT_RATES/STATIONS_ONLY_lon_lat_elev_gaugeID.nc'))
    }

    return(netcdf_file)
}

#'
#' Read all the gauges.
#'
get_all_gauges<-function(){

    netcdf_file = config_env$adjust_path_to_gdata_base_location(
        '/g/data/fj6/PTHA/AustPTHA_1/EVENT_RATES/STATIONS_ONLY_lon_lat_elev_gaugeID.nc')

    all_gauges = read_table_from_netcdf(netcdf_file)
    all_gauges$gaugeID = clean_gaugeID(all_gauges$gaugeID)

    return(all_gauges)
}

#'
#' Get lon/lat/depth/gaugeID from gauges in netcdf file
#' 
#' @param netcdf_file name of netdf tide gauge file
#' @param indices_of_subset integer vector with rows of the table to return. By
#' default return all rows
#' @param FASTREAD If true, check if we are trying to get a ptha18 file on the NCI
#' THREDDS server. If we are, then read instead from a "STATIONS_ONLY" file. This can
#' be quicker than reading the coordinates alone.
#' @return data.frame with lon, lat, elev, gaugeID
#' @export
get_netcdf_gauge_locations<-function(netcdf_file, indices_of_subset = NULL, FASTREAD=TRUE){

    if(FASTREAD){
        # Change the local copy of netcdf_file (in this function ONLY) to point to a file
        # that contains only station data.
        netcdf_file = get_file_with_gauges_only_if_on_NCI_THREDDS(netcdf_file)
    }
    fid = nc_open(netcdf_file, readunlim=FALSE)

    # Most files have the elevation in a variable 'elevation0',
    # but some have it in a variable 'elev'
    elev_var_name = 'elevation0'
    if(length(grep(elev_var_name, names(fid$var))) == 0){
        elev_var_name = 'elev'
        if(length(grep(elev_var_name, names(fid$var))) == 0){
            nc_close(fid)
            stop(paste0('No elevation variable found in file ', netcdf_file))
        }
    }

    if(is.null(indices_of_subset)){
        # Read everything
        gauge_ids = clean_gaugeID(ncvar_get(fid, 'gaugeID'))
        lon = ncvar_get(fid, 'lon')
        lat = ncvar_get(fid, 'lat')
        elev = ncvar_get(fid, elev_var_name)

    }else{

        # Read in one chunk -- faster over a remote filesystem
        contiguous_inds = min(indices_of_subset):max(indices_of_subset)
        contig_match = match(indices_of_subset, contiguous_inds)

        # In case it is faster to read everything -- seems to be over a remote connection
        gauge_ids = clean_gaugeID(ncvar_get(fid, 'gaugeID', start=contiguous_inds[1], count=length(contiguous_inds)))
        gauge_ids = gauge_ids[contig_match]
        lon = ncvar_get(fid, 'lon', start=contiguous_inds[1], count=length(contiguous_inds))
        lon = lon[contig_match]
        lat = ncvar_get(fid, 'lat', start=contiguous_inds[1], count=length(contiguous_inds))
        lat = lat[contig_match]
        elev = ncvar_get(fid, elev_var_name, start=contiguous_inds[1], count=length(contiguous_inds))
        elev = elev[contig_match]

    }

    nc_close(fid)

    gauge_data = data.frame(lon=lon, lat=lat, elev=elev, gaugeID = gauge_ids)

    return(gauge_data)
}

#' Extract the time vector from a netcdf file
#'
get_netcdf_gauge_output_times<-function(netcdf_file){

    fid = nc_open(netcdf_file, readunlim=FALSE)
    times = ncvar_get(fid, 'time')
    nc_close(fid)
    return(times)
}

#'
#' Find the indices of gauges in a netcdf file which are nearest to 
#' a matrix of lon/lat coordinates.
#'
#' @param netcdf_file name of netcdf tide gauge file
#' @param lonlat matrix with 2 columns containing lon/lat coordinates along rows
#' @return integer vector, containing the indices of points in the tide gauge
#' file which are nearest each point of lonlat
#' @export
get_netcdf_gauge_indices_near_points<-function(netcdf_file, lonlat){

    gauge_data = get_netcdf_gauge_locations(netcdf_file)
    lg = length(gauge_data[,1])
  
    inds = lonlat_nearest_neighbours(lonlat, cbind(gauge_data$lon, gauge_data$lat))

    return(inds)
}

#'
#' This code tests that netcdf gauge indices near points is working
#' 
.test_get_netcdf_gauge_indices_near_points<-function(netcdf_file){

    gauge_points = get_netcdf_gauge_locations(netcdf_file)

    test_indices = sample(1:length(gauge_points[,1]), size=10, replace=TRUE)

    derived_indices = get_netcdf_gauge_indices_near_points(netcdf_file, 
        gauge_points[test_indices,1:2])

    # It is possible that derived_indices != test_indices, if there are
    # repeated lon/lat points in the gauges
    #if(all(derived_indices == test_indices)){
    if(all(gauge_points[derived_indices,1:2] == gauge_points[test_indices, 1:2])){
        print('PASS')
    }else{
        print('FAIL')
    }

}

#' Find indices of gauges in netcdf file which are inside a polygon
#'
#' @param netcdf_file file with netcdf tide gauge data
#' @param region_poly SpatialPolygons object, inside which we want to know the
#' indices of all gauges. Alternatively it can be a 2 column matrix with x/y defining the polygon
#' @param return integer vector with indices of gauges inside the polygon
#' @export
get_netcdf_gauge_indices_in_polygon<-function(netcdf_file, region_poly){

    if(is.matrix(region_poly)){
        # If region_poly is a matrix, assume it is giving lon/lat coordinates
        # Then convert to a polygon below

        if(!((ncol(region_poly) == 2) & (nrow(region_poly) > 2))){
            print(region_poly)
            stop('If region_poly is a matrix, it must have 2 columns and more than 2 rows, in order to represent a polygon')
        }

        rp_dim = dim(region_poly)

        # Close the polygon if it is not closed
        if(!all(region_poly[1,] == region_poly[rp_dim[1],])){
            region_poly = rbind(region_poly, region_poly[rp_dim[1],])
        }

        region_poly = SpatialPolygons(list(Polygons(list(
            Polygon(region_poly, hole=FALSE)), ID=1)),
            proj4string=CRS("+init=epsg:4326"))
    }

    gauge_points = get_netcdf_gauge_locations(netcdf_file)

    coords = cbind(gauge_points$lon, gauge_points$lat)
    coords_sp = SpatialPoints(coords, 
        proj4string=CRS(proj4string(region_poly)))

    indices_of_subset = which(!is.na(over(coords_sp, region_poly)))
    if(length(indices_of_subset) == 0) stop('No coordinates in region_poly') 

    return(indices_of_subset)
}

#' Test of the above routine.
.test_get_netcdf_gauge_indices_in_polygon<-function(netcdf_file){

    gauge_points = get_netcdf_gauge_locations(netcdf_file)

    test_inds = c(1, length(gauge_points[,1]))

    test_poly = gBuffer(SpatialPoints(gauge_points[test_inds,1:2]), width=1.0e-03)
    proj4string(test_poly) = '+init=epsg:4326'
   
    inside_points = get_netcdf_gauge_indices_in_polygon(netcdf_file, test_poly)

    if(all(sort(test_inds) == sort(inside_points))){
        print('PASS')
    }else{
        print('FAIL')
    }
}

#' Find the index gauges in netcdf_file that match the given gauge_ID's
#'
#' @param netcdf_file netcdf filename
#' @param gauge_ID vector of numeric gauge IDs
#' @param FASTREAD if possible, read from a file that ONLY contains station data
#'
get_netcdf_gauge_index_matching_ID<-function(netcdf_file, gauge_ID, FASTREAD=TRUE){

    fid = nc_open(netcdf_file, readunlim=FALSE)

    if(FASTREAD){
        netcdf_stations_only = get_file_with_gauges_only_if_on_NCI_THREDDS(netcdf_file)
        fid2 = nc_open(netcdf_stations_only, readunlim=FALSE)
        point_ids = clean_gaugeID(ncvar_get(fid2, 'gaugeID'))
        nc_close(fid2)
    }else{
        point_ids = clean_gaugeID(ncvar_get(fid, 'gaugeID'))
    }

    closest_ID = match(gauge_ID, point_ids)

    if(any(is.na(closest_ID))){
        kk = which(is.na(closest_ID))
        print(gauge_ID[kk])
        nc_close(fid)
        stop('Could not find gauge index matching the above value (which is rounded to address netcdf issues)')
    }

    nc_close(fid)

    return(closest_ID)
}


.test_get_netcdf_gauge_index_matching_ID<-function(netcdf_file){

    netcdf_file = get_file_with_gauges_only_if_on_NCI_THREDDS(netcdf_file)

    desired_id = c(1.1, 10.1)
    desired_index = get_netcdf_gauge_index_matching_ID(netcdf_file, desired_id)

    # Check that the gaugeID value is 1.1 (to within precision of a float)
    fid = nc_open(netcdf_file)
    gaugeIDs = clean_gaugeID(ncvar_get(fid, 'gaugeID'))

    if(all(abs(gaugeIDs[desired_index] - desired_id) < 1.0e-08)){
        print('PASS')
    }else{
        print('FAIL')
    }
    nc_close(fid)

}


#' Convenience function to sort the netcdf tide gauge files in the same order as
#' the unit_source_statistics table. This is required for the unit-source
#' summation.
#'
#' @param unit_source_statistics data.frame with unit_source_statistics summary information,
#' typically derived as output from \code{rptha::discretized_source_summary_statistics}
#' @param netcdf_tide_gauge_files vector of netcdf tide gauge files. There should be at least one
#' file corresponding to each row of the unit_source_statistics. If there are more than one, then
#' only the first is used. There can also be extraneous files, which will be ignored.
#' @param character vector of netcdf_tide_gauge_files, ordered to correspond to
#' rows of unit_source_statistics
#' @export
sort_tide_gauge_files_by_unit_source_table<-function(
    unit_source_statistics,
    netcdf_tide_gauge_files){

    tide_gauge_rasters = sapply(netcdf_tide_gauge_files, 
        get_netcdf_attribute_initial_stage_raster)

    unit_source_to_tg = match(unit_source_statistics$initial_condition_file, 
        tide_gauge_rasters)

    if(any(is.na(unit_source_to_tg))){
        stop('Some initial condition files are not matching a corresponding tide gauge netcdf attribute')
    }

    return(netcdf_tide_gauge_files[unit_source_to_tg])
}

#' Read chunk of netcdf file
#'
#' Function for efficiently reading in a subset of the netcdf file 2D variables
#' (stage/uh/vh)
#'
#' @param fid netcdf file id
#' @param varname variable name to read
#' @param indices_of_subset desired indices for the second dimension (read all th first dimension)
#' @param chunk_size Read contiguous chunks of size 'chunk_size'. Over the internet it can be 
#' much faster to have this > 1, whereas on a disk it may be faster to be equal to one.
#' 
.chunked_read<-function(fid, varname, indices_of_subset, chunk_size=1){

    # Find the distance of every index from the first one
    di = indices_of_subset - min(indices_of_subset) 

    # Get indices in groups of 'chunk-size' from first offset
    # FIXME: This approach could be improved, by instead focussing on 
    # the 'jumps' between ordered indices_of_subset values.
    di_round = floor(di/chunk_size)
    unique_di_round = unique(di_round)

    counter = 0
    for(u in unique_di_round){
        # We will read all these gauges at once
        inds_to_get = which(di_round == u) 
        inds_range = range(indices_of_subset[inds_to_get])
        # Make a contiguous range covering the indices to read
        contiguous_inds = inds_range[1]:inds_range[2]

        # Find the indices we actually want in this contiguous range
        index_match = match(indices_of_subset[inds_to_get], contiguous_inds)

        # Read the data
        tmp = t(ncvar_get(fid, varname, 
            start=c(1,contiguous_inds[1]),
            count=c(-1, length(contiguous_inds))))

        if(counter == 0){
            values = matrix(NA, nrow=length(indices_of_subset), ncol=ncol(tmp))
        }
        counter = counter+1

        # Pack it in the output matrix
        values[inds_to_get,] = tmp[index_match,,drop=FALSE]
    }

    return(values)

}

#' Read flow time-series from the tide gauges netcdf file. 
#'
#' Optionally only read a subset of the flow time-series. Note that depending on 
#' the order of dimensions in the netcdf file, it may be more
#' time-efficient to read all the flow variables and then subset, or to just read
#' the subset
#'
#' @param netcdf_file filename
#' @param indices_of_subset if not NULL, then a vector of integer indices
#' corresponding to the gauges at which flow values are extracted
#' @param flow_and_attributes_only logical. If TRUE, then only read the flow
#' and the file attributes (other variables will be set to NULL). This can be
#' useful for efficiency
#' @param all_flow_variables logical. If TRUE, return an array
#' flow_time_series[i,j,k] with stage, uh, vh corresponding to k=1, 2, 3 resp. 
#' Otherwise return a matrix flow_time_series[i,j] containing just stage.
#' @return list containing a matrix of flow values, netcdf attributes,
#' gaugeIDs, and the value of indices_of_subset. FIXME: Consider returning
#' elevation as well
#' @export
#'
get_flow_time_series_SWALS<-function(netcdf_file, indices_of_subset=NULL,  
    flow_and_attributes_only=TRUE, all_flow_variables=TRUE){

    fid = nc_open(netcdf_file, readunlim=FALSE)
    # Determine if we are reading a remote file -- in that case file.exists will
    # return FALSE
    is_remote = (!file.exists(netcdf_file))
    # Get all 'global' attributes in the netcdf file. This allows sanity checks
    # on book_keeping, since the attributes include e.g. references to the unit
    # source initial conditions.
    run_atts = ncatt_get(fid, varid=0)

    if(!flow_and_attributes_only){
        netcdf_file2 = get_file_with_gauges_only_if_on_NCI_THREDDS(netcdf_file)
        fid2 = nc_open(netcdf_file)
        gauge_ids = clean_gaugeID(ncvar_get(fid2, 'gaugeID'))
        lat = ncvar_get(fid2, 'lat')
        lon = ncvar_get(fid2, 'lon')
        elev = ncvar_get(fid2, 'elevation0')
        nc_close(fid2)
    }else{
        # Don't bother reading other gauge information
        gauge_ids = NULL
        lat = NULL
        lon = NULL
        elev = NULL
    }

    # Get the names of the dimensions of stage -- this helps us determine
    # whether we can read the stage efficiently at single stations
    stage_dim1 = fid$var$stage$dim[[1]]$name
    stage_dim2 = fid$var$stage$dim[[2]]$name

    if(stage_dim1 == 'station'){
        # Read all stages, even if we only want a subset, because experience
        # suggests it is faster. This is because time is the 'slowly varying'
        # dimension, so we have to read through the entire file regardless.
        stages = ncvar_get(fid, 'stage')
        # Do subsetting if required
        if(!is.null(indices_of_subset)){
            stages = stages[indices_of_subset,,drop=FALSE]
            if(!flow_and_attributes_only){
                gauge_ids = gauge_ids[indices_of_subset]
                lon = lon[indices_of_subset]
                lat = lat[indices_of_subset]
                elev = elev[indices_of_subset]
            }
        }

        if(all_flow_variables){
            # Get x flux variable
            uhs = ncvar_get(fid, 'uh')
            if(!is.null(indices_of_subset)) uhs = uhs[indices_of_subset,,drop=FALSE]

            # Get y flux variable
            vhs = ncvar_get(fid, 'vh')
            if(!is.null(indices_of_subset)) vhs = vhs[indices_of_subset,,drop=FALSE]
        }

    }else{
        # In this case, the file has time varying quickly, which permits
        # efficient single-station access, so long as indices_of_subset is
        # continuous
        stopifnot(stage_dim1 == 'time')

        if(is.null(indices_of_subset)){
            # We want to read all stages
            read_all_stages = TRUE
        }else{
            # We only want a subset of the stages, so we might
            # be able to read efficiently
            read_all_stages = FALSE
        }

        if(read_all_stages){
            # Take the transpose of the stages for consistency with the case
            # when time is an unlimited dimension
            stages = t(ncvar_get(fid, 'stage'))
            if(!is.null(indices_of_subset)) stages = stages[indices_of_subset,,drop=FALSE]

            if(all_flow_variables){
                # Get x flux variable
                uhs = ncvar_get(fid, 'uh')
                if(!is.null(indices_of_subset)) uhs = uhs[indices_of_subset,,drop=FALSE]

                # Get y flux variable
                vhs = ncvar_get(fid, 'vh')
                if(!is.null(indices_of_subset)) vhs = vhs[indices_of_subset,,drop=FALSE]
            }
        }else{
            # Efficient reading is possible if indices_of_subset is
            # continuous (or close to). Take the transpose of the stages for
            # consistency with the case when time is an unlimited dimension

            chunk_size = ifelse(is_remote, 100, 1) 
            #print(paste0('chunk size: ', chunk_size))
            #print('stage...')
            stages = .chunked_read(fid, 'stage', indices_of_subset, chunk_size=chunk_size)

            if(all_flow_variables){
                #print('uhs...')
                uhs = .chunked_read(fid, 'uh', indices_of_subset, chunk_size=chunk_size)
                #print('vhs...')
                vhs = .chunked_read(fid, 'vh', indices_of_subset, chunk_size=chunk_size)
            }
                
        }

        # We still might have to subset the other variables
        if((!flow_and_attributes_only)&(!is.null(indices_of_subset))){
            gauge_ids = gauge_ids[indices_of_subset]
            lon = lon[indices_of_subset]
            lat = lat[indices_of_subset]
        }

    }
    nc_close(fid)

    if(all_flow_variables){
        flow_time_series = array(NA, dim=c(dim(stages), 3))
        flow_time_series[,,1] = stages
        flow_time_series[,,2] = uhs
        flow_time_series[,,3] = vhs
    }else{
        flow_time_series = array(NA, dim=c(dim(stages), 1))
        flow_time_series[,,1] = stages
    }

    output = list(attr = run_atts, flow_time_series = flow_time_series, gauge_ids = gauge_ids, 
        indices_of_subset=indices_of_subset, lon=lon, lat=lat, elev=elev)
    return(output)
}

#' Combine tsunami unit sources into tsunami events
#'
#' This function combines data describing earthquake events with data describing
#' unit source geometries, and files containing flow time-series for each 
#' unit source. It can be used to produce flow time-series for earthquake events
#' (assuming linearity), or derive summary statistics from them. The user
#' provides the function to read the flow time-series (this will vary depending
#' on the tsunami propagation solver used). 
#' 
#' @param earthquake_events data.frame with earthquake events, e.g. from the 
#' output of \code{get_all_earthquake_events}
#' @param unit_source_statistics data.frame with unit source summary stats, e.g.
#' from the output of \code{discretized_source_summary_statistics}. It must
#' have a column named 'event_index_string' giving the unit-source indices in the
#' event, and a column 'slip' giving the earthquake slip.
#' @param unit_source_flow_files vector of filenames, corresponding to
#' rows of unit_source_statistics, containing flow time-series for each unit
#' source
#' @param get_flow_time_series_function function to read the flow_time_series
#' data. This function is provided by the user, and is expected to vary depending
#' on which flow solver is used. The function must take as input a single
#' unit_source_flow_file, and an optional indices_of_subset vector which
#' contains the indices of flow_time_series to extract, or NULL to extract all
#' indices. This function must returns a list containing an element named
#' "flow_time_series", which is EITHER an array flow_time_series[i,j,k] containing stage
#' time-series, with i indexing the gauge, and j the time stage time slice, and k=1, OR
#' an array including variables other than just stage, such as uh, vh. In this case k = 1,2,...
#' (e.g. flow_time_series[i,j,k] where i indexes the station, j the time, and k
#' the variable, typically k=1 is stage, k=2 is uh, and k=3 is vh). The list
#' that the function returns can also contain entries with other names, but they
#' are ignored by this function
#' @param indices_of_subset If not null, the value of indices_of_subset that is
#' passed to get_flow_time_series_function
#' @param verbose logical. If TRUE, print information on progress
#' @param summary_function function which takes the flow data for a single 
#' earthquake event as input, and returns 'anything you like' as output. This is applied 
#' to the flow data for each event before it is output. If NULL, the full flow
#' time-series is provided. This function could be used to compute e.g.
#' the maxima and period of the stage time-series at each gauge, while avoiding
#' having to store the full waveforms for all events in memory.
#' @param msl mean sea level for the linear shallow water solver. All gauges above msl
#' are 'dry' for the full simulation
#' @param ... further arguments passed to get_flow_time_series_function
#' @export
#'
make_tsunami_event_from_unit_sources<-function(
    earthquake_events, 
    unit_source_statistics, 
    unit_source_flow_files,
    get_flow_time_series_function = get_flow_time_series_SWALS,  
    indices_of_subset=NULL, 
    verbose=FALSE,
    summary_function=NULL,
    msl=0.0,
    ...){

    if(all(c('slip', 'event_slip_string') %in% names(earthquake_events))){
        msg = paste0('earthquake_events cannot have both a column named "slip" ', 
            'and a column named "event_slip_string", since the presence of ', 
            'one or the other is used to distinguish stochastic slip from uniform slip')
        stop(msg)
    }

    if(verbose) print('Finding which unit sources we need ...')

    num_eq = length(earthquake_events[,1])
    events_data = vector(mode='list', length=num_eq)

    # Figure out which unit sources we need flow_time_series for 
    required_unit_sources = vector(mode='list', length=num_eq)
    for(i in 1:num_eq){
        required_unit_sources[[i]] = 
            get_unit_source_indices_in_event(earthquake_events[i,])
    }

    union_required_unit_sources = unique(unlist(required_unit_sources))

    # Double check that the unit_source_statistics are correctly numbered, and
    # ordered
    stopifnot(all(unit_source_statistics$subfault_number == 
        1:length(unit_source_statistics$subfault_number)))

    if(verbose) print('Reading required unit sources ...')

    # Read all the unit source tsunami data that we require
    # We make a list that is long enough to hold all unit source tsunami
    # results , but only read in the ones that we need.
    flow_data = vector(mode='list', 
        length=length(unit_source_flow_files))
    for(i in union_required_unit_sources){

        netcdf_file = unit_source_flow_files[i]

        flow_data[[i]] = get_flow_time_series_function(
            netcdf_file, 
            indices_of_subset = indices_of_subset,
            ...)

        if(msl != 0){
            # Subtract MSL from stage, since (stage-msl) is linear 
            flow_data[[i]]$flow_time_series[,,1] = flow_data[[i]]$flow_time_series[,,1] - msl
        }

        names(flow_data)[i] = netcdf_file

        #gc()
    }

    if(verbose) print('Summing unit sources ...') 

    # Make a matrix in which we sum the flow_time_series for each event
    first_uss = which(!is.na(names(flow_data)))[1]
    template_flow_data = flow_data[[first_uss]]$flow_time_series * 0

    if('slip' %in% names(earthquake_events)){
        uniform_slip = TRUE
    }else{
        if(!('event_slip_string' %in% names(earthquake_events))){
            stop('earthquake_events must either have a column named "slip", OR one named "event_slip_string"')
        }
        uniform_slip = FALSE
    }

    # Do the sum
    for(i in 1:num_eq){

        if(verbose) print(paste0('    event ', i))
        earthquake_event = earthquake_events[i,]
        event_unit_sources = required_unit_sources[[i]]

        template_flow_data = template_flow_data * 0
        #
        # NOTE: We need to be sure that R is not treating "template_flow_data"
        # as a reference to some other variable, since below we update in-place
        # with "axpy_local".
        #
        # For example, in R, the following 'axpy_local' call will modify both
        # 'x' and 'z', because R will not have made a copy for z [instead, x
        # and z point to the same memory]:
        #     x = runif(10)
        #     z = x
        #     y = runif(10)
        #     axpy_local(x, 1.0, y)
        # On the other hand, we can force a copy by creating z as:
        #     z = x*1.0   # The multiplication forces a copy
        # and then the problem does not arise
        #

        if(uniform_slip){
            # Sum the unit sources [each with 1m slip]
            for(j in event_unit_sources){
                ## This is slower than BLAS
                #template_flow_data = template_flow_data + 
                #    flow_data[[j]]$flow_time_series

                ## Does the same as commented out above, faster with BLAS
                axpy_local(template_flow_data, 1.0, flow_data[[j]]$flow_time_series)
            }

            template_flow_data = template_flow_data * earthquake_event$slip

        }else{
            # Stochastic slip case

            # read the slip vector in it's funny character format.
            slip_vector = scan(text=gsub('_', ' ', earthquake_event$event_slip_string), quiet=TRUE)
            stopifnot(length(slip_vector) == length(event_unit_sources))

            # Sum the non-uniform slip
            for(j in 1:length(event_unit_sources)){
                ## Slower than BLAS
                #template_flow_data = template_flow_data +
                #    flow_data[[event_unit_sources[j]]]$flow_time_series * slip_vector[j]

                # Does the same as commented out above, with BLAS
                ej = event_unit_sources[[j]]
                axpy_local(template_flow_data, slip_vector[j], flow_data[[ej]]$flow_time_series)
            }

        }

        # Add MSL to stage [recall it was subtracted earlier, since (stage-msl) is linear]
        if(msl != 0) template_flow_data[,,1] = template_flow_data[,,1] + msl

        # Apply the summary function if required
        if(is.null(summary_function)){
            events_data[[i]] = template_flow_data
        }else{
            events_data[[i]] = summary_function(template_flow_data)
        }

    }

    rm(flow_data, template_flow_data)
    gc()

    return(events_data)
}


#
#' Generic routine to run the test codes on a netcdf file
#'
#' Helps to detect connection issues as well as errors in code logic
#' 
#' Note the test_all.R code also runs a test on actual waveform summation
#' (whereas the tests below only check gauge-finding routines -- which is limited).
test_sum_tsunami_unit_sources<-function(netcdf_file){

    .test_get_netcdf_gauge_indices_near_points(netcdf_file)
    .test_get_netcdf_gauge_indices_in_polygon(netcdf_file)
    .test_get_netcdf_gauge_index_matching_ID(netcdf_file)

}
