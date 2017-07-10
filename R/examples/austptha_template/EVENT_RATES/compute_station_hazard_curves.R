library(rptha)
library(parallel)
#
# INPUTS
#

# Files with uniform slip max_stage for every point, and also event rates
all_source_uniform_slip_tsunami = Sys.glob(
    '../SOURCE_ZONES/*/TSUNAMI_EVENTS/all_uniform_slip_earthquake_events_tsunami_*.nc')
# Files with uniform slip max_stage for every point, and also event rates
all_source_stochastic_slip_tsunami = Sys.glob(
    '../SOURCE_ZONES/*/TSUNAMI_EVENTS/all_stochastic_slip_earthquake_events_tsunami_*.nc')

# Apply hazard curve computation / data extraction to chunks of data this big
point_chunk_size = 100

# Stages at which we compute the rate, for every point, for every source-zone
stage_seq_max = 20 # m
stage_seq_min = 0.02 # m
stage_seq_len = 100 # Number of points
stage_seq = exp(seq(log(stage_seq_min), log(stage_seq_max), len=stage_seq_len))

# Output file name


#
# End inputs
#


#'
#'  Get lon/lat/elevation for all gauges
#' 
read_lon_lat_elev<-function(nc_file){

    fid = nc_open(nc_file, readunlim=FALSE)
   
    lon = ncvar_get(fid, 'lon') 
    lat = ncvar_get(fid, 'lat') 
    elev = ncvar_get(fid, 'elev')
    gaugeID = ncvar_get(fid, 'gaugeID')

    nc_close(fid)

    output = data.frame(lon=lon, lat=lat, elev=elev, gaugeID = gaugeID)

    return(output)
}

#'
#' Get stage exceedance rates
#' 
#' @param tsunami_file netcdf file with tsunami max_stage / rates, for every gauge and event
#' @param gauge_points data.frame with gauge coordinates and elevation, to check that point
#'     ordering in tsunami_file is as expected 
#' @param point_chunk_size Number of gauges for which we read max_stage simultaneously
#' @param stage_seq stages at which we return the integrated exceedance rate
source_zone_stage_exceedance_rates<-function(tsunami_file, gauge_points, point_chunk_size, stage_seq){

    # Check that gauges are ordered as expected
    lgp = length(gauge_points[,1])
    local_gauge_points = read_lon_lat_elev(tsunami_file)
    stopifnot(all(local_gauge_points == gauge_points))

    # Store the rates at every stage in stage_seq
    output_rates = matrix( NA, nrow=stage_seq_len, ncol=lgp )
    output_rates_lower_ci = matrix( NA, nrow=stage_seq_len, ncol=lgp )
    output_rates_upper_ci = matrix( NA, nrow=stage_seq_len, ncol=lgp )

    # Process gauges in chunks to reduce memory usage
    point_chunks = splitIndices(lgp, ceiling(lgp/point_chunk_size))

    # File we read from
    fid = nc_open(tsunami_file, readunlim=FALSE) 

    # Get rate for every event
    event_rate = ncvar_get(fid, 'event_rate_annual')
    event_rate_lower = ncvar_get(fid, 'event_rate_annual_lower_ci')
    event_rate_upper = ncvar_get(fid, 'event_rate_annual_upper_ci')

    # Do the calculation for all points, in chunks
    for(j in 1:length(point_chunks)){

        #print(paste0(c(j, '/', length(point_chunks))))

        gauge_indices = point_chunks[[j]]
        # Ensure indices are contiguous
        stopifnot(all(range(diff(gauge_indices)) == 1))

        # Get max stage for every event, for this chunk of gauges
        peak_stages = ncvar_get(fid, 'max_stage', 
            start=c(1, gauge_indices[1]), 
            count=c(-1, length(gauge_indices)))

        for(k in 1:ncol(peak_stages)){

            if(all(is.na(peak_stages[,k]))) next
           
            # Sort the stages in decreasing order 
            events_sort = sort(peak_stages[,k], index.return=TRUE, decreasing=TRUE)
            gauge_no = gauge_indices[k]

            # Value that no stage is larger than for this gauge
            stage_upper_bound = events_sort$x[1] + 1.0e-03
            # Value no stage is larger than for ANY gauge!
            global_stage_upper_bound = 99999999

            # Sort the rates like the stages, and do cumulative sum to get an
            # exceedance rate
            sorted_stages = c(global_stage_upper_bound, stage_upper_bound, events_sort$x)
            sorted_cumulative_rate = cumsum(c(0, 0, event_rate[events_sort$ix]))
            sorted_cumulative_rate_lower = cumsum(c(0, 0, event_rate_lower[events_sort$ix]))
            sorted_cumulative_rate_upper = cumsum(c(0, 0, event_rate_upper[events_sort$ix]))

            # Rates
            rate_fun = approx(sorted_stages, sorted_cumulative_rate, xout=stage_seq)
            output_rates[,gauge_no] = rate_fun$y

            # Upper estimate of rates
            rate_fun = approx(sorted_stages, sorted_cumulative_rate_upper, xout=stage_seq)
            output_rates_upper_ci[,gauge_no] = rate_fun$y
                        
            # Lower estimate of rates
            rate_fun = approx(sorted_stages, sorted_cumulative_rate_lower, xout=stage_seq)
            output_rates_lower_ci[,gauge_no] = rate_fun$y
        }
    }

    output = list(
        rates = output_rates, 
        rates_upper_ci = output_rates_upper_ci, 
        rates_lower_ci = output_rates_lower_ci)

    return(output)
}

#' Take care of saving outputs to netcdf file
#'
#'
create_rate_netcdf_file<-function(
    source_name, 
    gauge_points, 
    stage_seq, 
    uniform_slip_rates, 
    stochastic_slip_rates,
    uniform_slip_tsunami_file, 
    stochastic_slip_tsunami_file){

    # Dimension for rate curve
    dim_stage_seq = ncdim_create('stage', 'm', vals=stage_seq, unlim=FALSE,
        longname='stages corresponding to tsunami wave height exceedance rates')

    # Dimension for gauges
    dim_station = ncdim_create('station', '', vals=1:length(gauge_points[,1]), 
        unlim=TRUE,
        longname='integer index corresponding to the gauge location')

    # Variables for gauge locations
    gauge_lon_v = ncvar_def(name='lon', units='degrees_east', dim=list(dim_station), 
        missval=NA, longname='station_longitude', prec='float')
    gauge_lat_v = ncvar_def(name='lat', units='degrees_north', dim=list(dim_station), 
        missval=NA, longname='station_latitude', prec='float')
    gauge_elev_v = ncvar_def(name='elev', units='m', dim=list(dim_station), 
        missval=NA, longname='station_ground_elevation_above_mean_sea_level', 
        prec='float')
    gauge_id_v = ncvar_def(name='gaugeID', units='', dim=list(dim_station), 
        missval=NA, longname='real_ID_for_each_station', prec='float')

    all_nc_var = list(gauge_lon_v, gauge_lat_v, gauge_elev_v, gauge_id_v)

    # Variables for rates, uniform slip
    uniform_rate_v = ncvar_def(
        name='uniform_slip_rate', units='events per year',
        dim=list(dim_stage_seq, dim_station), 
        longname = 'exceedance rate of peak stage for uniform slip events')

    uniform_rate_upper_v = ncvar_def(
        name='uniform_slip_rate_upper_ci', units='events per year',
        dim=list(dim_stage_seq, dim_station), 
        longname = 'exceedance rate (upper credible interval) of peak stage for uniform slip events')

    uniform_rate_lower_v = ncvar_def(
        name='uniform_slip_rate_lower_ci', units='events per year',
        dim=list(dim_stage_seq, dim_station), 
        longname = 'exceedance rate (lower credible interval) of peak stage for uniform slip events')

    all_nc_var = c(all_nc_var, list(uniform_rate_v, uniform_rate_upper_v, uniform_rate_lower_v))

    # Variables for rates, stochastic slip
    stochastic_rate_v = ncvar_def(
        name='stochastic_slip_rate', units='events per year',
        dim=list(dim_stage_seq, dim_station), 
        longname = 'exceedance rate of peak stage for stochastic slip events')

    stochastic_rate_upper_v = ncvar_def(
        name='stochastic_slip_rate_upper_ci', units='events per year',
        dim=list(dim_stage_seq, dim_station), 
        longname = 'exceedance rate (upper credible interval) of peak stage for stochastic slip events')

    stochastic_rate_lower_v = ncvar_def(
        name='stochastic_slip_rate_lower_ci', units='events per year',
        dim=list(dim_stage_seq, dim_station), 
        longname = 'exceedance rate (lower credible interval) of peak stage for stochastic slip events')

    all_nc_var = c(all_nc_var, list(stochastic_rate_v, stochastic_rate_upper_v, stochastic_rate_lower_v))

    # Make name for output file
    sourcename_dot_nc = paste0(source_name, '.nc')
    output_file_name = paste0(
        dirname(uniform_slip_tsunami_file),
        '/', 
        'tsunami_stage_exceedance_rates_', sourcename_dot_nc)

    # Create output file
    output_fid = nc_create(output_file_name, vars=all_nc_var)

    # Attributes
    ncatt_put(output_fid, varid=0, attname = 'uniform_slip_tsunami_event_file',
        attval=normalizePath(uniform_slip_tsunami_file), prec='text')
    ncatt_put(output_fid, varid=0, attname = 'stochastic_slip_tsunami_event_file',
        attval=normalizePath(stochastic_slip_tsunami_file), prec='text')
    ncatt_put(output_fid, varid=0, attname='source_zone_name',
        attval=source_name, prec='text')
    ncatt_put(output_fid, varid=0, attname='parent_script_name',
        attval=parent_script_name(), prec='text')

    # Values of gauges
    ncvar_put(output_fid, gauge_lon_v, gauge_points$lon) 
    ncvar_put(output_fid, gauge_lat_v, gauge_points$lat) 
    ncvar_put(output_fid, gauge_elev_v, gauge_points$elev) 
    ncvar_put(output_fid, gauge_id_v, gauge_points$gaugeID) 

    # Values of uniform rates
    ncvar_put(output_fid, uniform_rate_v, uniform_slip_rates$rates)
    ncvar_put(output_fid, uniform_rate_upper_v, uniform_slip_rates$rates_upper_ci)
    ncvar_put(output_fid, uniform_rate_lower_v, uniform_slip_rates$rates_lower_ci)

    # Values of stochastic rates
    ncvar_put(output_fid, stochastic_rate_v, stochastic_slip_rates$rates)
    ncvar_put(output_fid, stochastic_rate_upper_v, stochastic_slip_rates$rates_upper_ci)
    ncvar_put(output_fid, stochastic_rate_lower_v, stochastic_slip_rates$rates_lower_ci)

    nc_close(output_fid)

    return(invisible(output_file_name))
}

stop()
# Get point info
gauge_points = read_lon_lat_elev(all_source_uniform_slip_tsunami[1])

# names of sources
source_names = basename(dirname(dirname(all_source_uniform_slip_tsunami)))
# check that uniform/stochastic files are ordered the same
stopifnot(all(source_names == basename(dirname(dirname(all_source_stochastic_slip_tsunami)))  ) )

for(i in 1:length(source_names)){

    # source information
    source_name = source_names[i]
    uniform_slip_tsunami_file = all_source_uniform_slip_tsunami[i]
    stochastic_slip_tsunami_file = all_source_stochastic_slip_tsunami[i]

    # Get uniform slip outputs
    uniform_slip_rates = source_zone_stage_exceedance_rates(
        uniform_slip_tsunami_file, gauge_points, 
        point_chunk_size=100, stage_seq=stage_seq)

    # Get stochastic slip outputs
    stochastic_slip_rates = source_zone_stage_exceedance_rates(
        stochastic_slip_tsunami_file, gauge_points, 
        point_chunk_size=100, stage_seq=stage_seq)

    create_rate_netcdf_file(
        source_name, gauge_points, stage_seq, 
        uniform_slip_rates, stochastic_slip_rates,
        uniform_slip_tsunami_file, stochastic_slip_tsunami_file)    

}

