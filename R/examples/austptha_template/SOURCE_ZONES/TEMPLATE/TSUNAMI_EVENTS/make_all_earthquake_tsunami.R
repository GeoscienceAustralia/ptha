library(rptha)
library(parallel)
source('sum_tsunami_unit_sources.R')

#
# Input parameters
#

# Memory we will allow in unit-sources, in MB. 
# It's a rough estimator -- should be a fraction of 1-nodes memory [e.g. using
# 10 GB on a 32GB machine seemed to work ok]
memory_for_unit_sources = 10 * 1024 

# Number of cores in parallel [shared memory only]
mc_cores = 16 

# Run parameters
msl = 0.0 # Pre-tsunami initial condition for linear model. Gauges with elev>msl are inactive
lat_range = c(-72 + 2/60, 65 - 2/60) # Latitude range where gauges should be taken
#lon_range = c(-40, 320) # Longitude range where gauges can be taken -- currently unused

# Record the wave arrival time as the time at which abs(stage) exceeds this
# threshold
stage_threshold_for_arrival_time = 0.001 

# Files to read
source_zone_name = basename(dirname(getwd()))
earthquake_events_file = paste0('all_uniform_slip_earthquake_events_', source_zone_name, '.nc')
unit_source_statistics_file = paste0('unit_source_statistics_', source_zone_name, '.nc')

#
# End input
#

# Read unit source statistics & uniform slip earthquake events
unit_source_statistics = read_table_from_netcdf(unit_source_statistics_file)
all_eq_events = read_table_from_netcdf(earthquake_events_file)

# Read gauges
gauge_locations = get_netcdf_gauge_locations(unit_source_statistics$tide_gauge_file[1])
# Read times at which gauge values are stored -- always constant thanks to the model setup
gauge_times = get_netcdf_gauge_output_times(unit_source_statistics$tide_gauge_file[1])

#
# Partition the work among processes
#
nevents = length(all_eq_events[,1])
ngauges = length(gauge_locations[,1])
nunit_sources = length(unit_source_statistics[,1])
real_bytes = 8  # size of an R real
nvar = 3 # [stage, uh, vh]
full_unit_sources_mem = real_bytes * nvar * length(gauge_times) * ngauges * 
    nunit_sources /(1024 * 1024)
# Make the station_chunk_size so that when spread over all cores, it 
# takes up memory = memory_for_unit_sources.
# My 'model' of how this code uses memory is very approximate and
# under-estimates, so you should set 'memory_for_unit_sources' to be say 1/3 of
# the available memory.
station_chunk_size = floor( length(gauge_locations[,1]) * 
    min(1, memory_for_unit_sources/(full_unit_sources_mem*mc_cores)) )
gauge_chunks_list = splitIndices(ngauges, ceiling(ngauges/station_chunk_size))

#
# Given an array with [gauge, time-slice, [stg,uh,vh]], extract flow summary
# statistics for each gauge. Pass this to the function that makes the tsunami
# events, so we don't have to store full time-series [too much memory for
# many gauges and events]
#
local_summary_function<-function(flow_data){

    times = gauge_times
    mean_dt = mean(diff(times))

    output = matrix(NA, nrow=dim(flow_data)[1], ncol=5)

    for(i in 1:dim(flow_data)[1]){
        stages = flow_data[i,,1]
        rst = range(stages)
        # Peak stage
        output[i,1] = max(rst)
        # Reference period
        # FIXME: This period is approximate in a number of ways -- consider revising
        #        e.g. start only after the wave exceeds some threshold?
        #             consider only the most significant wave?
        output[i,2] = rptha::zero_crossing_period(stages, dt=mean_dt)
        # FIXME: Peak_to_trough -- note it might be better to have a 'wave-based' view
        output[i,3] = diff(rst)
        # Arrival time
        suppressWarnings({kk = min(which(abs(stages) > stage_threshold_for_arrival_time))})
        if(is.finite(kk)){
            output[i,4] = times[kk]
        }

        # Initial stage
        output[i,5] = stages[1]

    }

    return(output)

}

#
# Loop over chunks of gauges, computing the statistics as we go
#

#modelled_flow_store = vector(mode='list', length=length(gauge_chunks_list))
#for(i in 1:length(gauge_chunks_list)){
#
#    gcl = gauge_chunks_list[[i]]
#
#    modelled_flow_store[[i]] = make_tsunami_event_from_unit_sources(
#        earthquake_events = all_eq_events,
#        unit_source_statistics = unit_source_statistics,
#        unit_source_flow_files = unit_source_statistics$tide_gauge_file,
#        indices_of_subset = gcl,
#        summary_function = local_summary_function,
#        msl=msl)
#}


# Parallel version of the above loop
parfun<-function(gcl){
    make_tsunami_event_from_unit_sources(
        earthquake_events = all_eq_events,
        unit_source_statistics = unit_source_statistics,
        unit_source_flow_files = unit_source_statistics$tide_gauge_file,
        indices_of_subset = gcl,
        summary_function = local_summary_function,
        msl = msl)
}
modelled_flow_store = mclapply(gauge_chunks_list, parfun, mc.cores=mc_cores)

#
# Pack outputs
#

# Matrix for max stage, with as many rows as events, and as many columns as gauges.
gauge_event_max_stage = matrix(-9999, nrow=nevents, ncol=ngauges)
# Period
gauge_event_reference_period = matrix(-9999, nrow=nevents, ncol=ngauges)
# Peak-to-trough
gauge_event_peak_to_trough = matrix(-9999, nrow=nevents, ncol=ngauges)
# Arrival time
gauge_event_arrival_time = matrix(-9999, nrow=nevents, ncol=ngauges)
# Initial stage
gauge_event_initial_stage = matrix(-9999, nrow=nevents, ncol=ngauges)

# Pack the statistics into the output arrays
for(i in 1:length(gauge_chunks_list)){
    gcl = gauge_chunks_list[[i]]

    # Vector with 1 at 'valid gauges', 0 elsewhere
    valid_gauges = (
        gauge_locations$lat[gcl] < lat_range[2] &
        gauge_locations$lat[gcl] > lat_range[1] &
        gauge_locations$elev[gcl] < msl)

    # Set 'invalid gauges' values to zero
    for(j in 1:nevents){
        gauge_event_max_stage[j, gcl ] = 
            modelled_flow_store[[i]][[j]][,1] * (c(NA, 1)[valid_gauges+1])
        gauge_event_reference_period[j, gcl ] = 
            modelled_flow_store[[i]][[j]][,2] * (c(NA, 1)[valid_gauges+1])
        gauge_event_peak_to_trough[j, gcl ] = 
            modelled_flow_store[[i]][[j]][,3] * (c(NA, 1)[valid_gauges+1])
        gauge_event_arrival_time[j, gcl ] = 
            modelled_flow_store[[i]][[j]][,4] * (c(NA, 1)[valid_gauges+1])
        gauge_event_initial_stage[j, gcl ] = 
            modelled_flow_store[[i]][[j]][,5] * (c(NA, 1)[valid_gauges+1])
    }
}



##############################################################################
#
# Write all the relevant information to netcdf
#
##############################################################################

write_all_source_zone_tsunami_statistics_to_netcdf(
    source_zone_name,
    all_eq_events,
    all_eq_events_file,
    unit_source_statistics,
    unit_source_statistics_file,
    gauge_locations,
    gauge_event_max_stage,
    gauge_event_reference_period,
    gauge_event_peak_to_trough,
    gauge_event_arrival_time,
    gauge_event_initial_stage,
    stage_threshold_for_arrival_time){

    library(ncdf4)

    # Station dimension -- make this one unlimited for fast 'single-station'
    # access to the data
    dim_station = ncdim_def(name='station', units='', 
        vals=1:length(gauge_locations[,1]), unlim=TRUE,
        longname='integer index corresponding to the gauge locations')

    # Event dimension
    dim_event = ncdim_def(name='event', units='', 
        vals=1:length(all_eq_events[,1]), unlim=FALSE,
        longname='integer index corresponding to the tsunami event')

    # Unit source dimension
    dim_unit_sources = ncdim_def(name='unitsource', units='', 
        vals=1:length(unit_source_statistics[,1]), unlim=FALSE,
        longname='integer index corresponding to the unit_sources')

    # Figure out a dimension size for a string
    charlen_sourcename = max(nchar(all_eq_events$sourcename))
    charlen_event_index_string = max(nchar(all_eq_events$event_index_string))
    charlen_initial_cond_string = max(nchar(unit_source_statistics$initial_condition_file))
    charlen_tide_gauge_string = max(nchar(unit_source_statistics$tide_gauge_file))

    dim_char_size = max(c(charlen_sourcename, charlen_event_index_string, 
        charlen_initial_cond_string, charlen_tide_gauge_string))

    dim_nchar = ncdim_def(name='max_nchar', units='', vals=1:dim_char_size, 
        unlim=FALSE,
        longname='integer index corresponding to the maximum number of string characters')


    #
    # Create netcdf variable for gauges
    #
    gauge_lon_v = ncvar_def(name='lon', units='degrees_east', dim=list(dim_station), 
        missval=NA, longname='station_longitude', prec='float')
    gauge_lat_v = ncvar_def(name='lat', units='degrees_north', dim=list(dim_station), 
        missval=NA, longname='station_latitude', prec='float')
    gauge_elev_v = ncvar_def(name='elev', units='m', dim=list(dim_station), 
        missval=NA, longname='station_ground_elevation_above_mean_sea_level', 
        prec='float')
    gauge_id_v = ncvar_def(name='gaugeID', units='', dim=list(dim_station), 
        missval=NA, longname='real_ID_for_each_station', prec='float')

    # Keep a list of variables for when we make the ncdf
    all_nc_var = list(gauge_lon_v, gauge_lat_v, gauge_elev_v, gauge_id_v)

    #
    # Create a netcdf variable for the gauge summary statistics
    #
    gauge_event_max_stage_v = ncvar_def(name='max_stage', units='m', 
        dim=list(dim_event, dim_station), missval=NA, 
        longname='maximum_stage_at_gauge_during_event',
        prec='float')
    gauge_event_reference_period_v = ncvar_def(name='period', units='s', 
        dim=list(dim_event, dim_station), missval=NA, 
        longname='zero_crossing_wave_period_of_gauge_time_series',
        prec='float')
    gauge_event_peak_to_trough_v = ncvar_def(name='stage_range', units='m', 
        dim=list(dim_event, dim_station), missval=NA, 
        longname='difference_between_the_maximum_stage_and_minimum_stage',
        prec='float')
    # Make a long name for the arrival time
    arrival_time_longname = paste0(
        'initial_time_at_which_stage_perturbation_exceeds_', 
        as.character(stage_threshold_for_arrival_time),
        '_meters_in_absolute_value')
    gauge_event_arrival_time_v = ncvar_def(name='arrival_time', units='s', 
        dim=list(dim_event, dim_station), missval=NA, 
        longname=arrival_time_longname,
        prec='float')
    gauge_event_initial_stage_v = ncvar_def(name='initial_stage', units='m', 
        dim=list(dim_event, dim_station), missval=NA, 
        longname='initial_stage_at_gauge',
        prec='float')

    # Keep a list of variables for when we make the netcdf 
    all_nc_var = c(all_nc_var, list(gauge_event_max_stage_v, 
        gauge_event_reference_period_v, gauge_event_peak_to_trough_v,
        gauge_event_arrival_time_v, gauge_event_initial_stage_v))

    #
    # Create a netcdf variable for the event summary statistics
    #
    event_area_v = ncvar_def(name='event_area', units='km^2', dim=list(dim_event),
        missval=NA, longname='rupture_area_of_earthquake_event', prec='float')
    event_mean_length_v = ncvar_def(name='event_mean_length', units='km', dim=list(dim_event),
        missval=NA, longname='mean_rupture_length_of_earthquake_event', prec='float')
    event_mean_width_v = ncvar_def(name='event_mean_width', units='km', dim=list(dim_event),
        missval=NA, longname='mean_rupture_width_of_earthquake_event', prec='float')
    event_slip_v = ncvar_def(name='event_slip', units='m', dim=list(dim_event),
        missval=NA, longname='slip_of_earthquake_event', prec='float')
    event_Mw_v = ncvar_def(name='event_Mw', units='', dim=list(dim_event),
        missval=NA, longname='moment_magnitude_of_earthquake_event', prec='float')
    event_mean_depth_v = ncvar_def(name='event_mean_depth', units='km', dim=list(dim_event),
        missval=NA, longname='mean_rupture_depth_of_earthquake_event', prec='float') 
    event_max_depth_v = ncvar_def(name='event_max_depth', units='km', dim=list(dim_event),
        missval=NA, longname='max_rupture_depth_of_earthquake_event', prec='float') 
    event_index_string_v = ncvar_def(name='event_index_string', units='', dim=list(dim_event, dim_nchar),
        missval=NA, longname='indices_of_unit_sources_included_in_earthquake_event', prec='char')
    event_sourcename_v = ncvar_def(name='event_sourcename', units='', dim=list(dim_event, dim_nchar),
        missval=NA, longname='source_zone_name_of_earthquake_event', prec='char')
    event_rate_annual_v = ncvar_def(name='event_rate_annual', units='events per year', dim=list(dim_event),
        missval=NA, longname='mean_rate_of_earthquake_event', prec='float') 
    event_rate_annual_upper_ci_v = ncvar_def(name='event_rate_annual_upper_ci', 
        units='events per year', dim=list(dim_event),
        missval=NA, longname='upper_credible_interval_for_rate_of_earthquake_event', 
        prec='float') 
    event_rate_annual_lower_ci_v = ncvar_def(name='event_rate_annual_lower_ci', 
        units='events per year', dim=list(dim_event),
        missval=NA, longname='lower_credible_interval_for_rate_of_earthquake_event', 
        prec='float')

    # Append to list of variables, for making netcdf
    all_nc_var = c(all_nc_var, list(event_area_v, event_mean_length_v, 
        event_mean_width_v, event_slip_v, event_Mw_v, event_mean_depth_v, 
        event_max_depth_v, event_index_string_v, event_sourcename_v, 
        event_rate_annual_v, event_rate_annual_upper_ci_v, 
        event_rate_annual_lower_ci_v) )

    #
    # Make variables for unit-source statistics
    #
    us_lon_c_v = ncvar_def(name='us_lon_c', units='degrees_east', 
        dim=list(dim_unit_sources), missval=NA, 
        longname='longitude_of_unit_source_centroid', prec='float')
    us_lat_c_v = ncvar_def(name='us_lat_c', units='degrees_north', 
        dim=list(dim_unit_sources), missval=NA, 
        longname='latitude_of_unit_source_centroid', prec='float')
    us_depth_v = ncvar_def(name='us_depth', units='km', 
        dim=list(dim_unit_sources), missval=NA, 
        longname='mean_depth_of_unit_source', prec='float')
    us_strike_v = ncvar_def(name='us_strike', units='degrees', 
        dim=list(dim_unit_sources), missval=NA, 
        longname='strike_of_unit_source_in_degrees_clockwise_from_north', 
        prec='float')
    us_dip_v = ncvar_def(name='us_dip', units='degrees', 
        dim=list(dim_unit_sources), missval=NA, 
        longname='dip_of_unit_source_in_degrees', 
        prec='float')
    us_rake_v = ncvar_def(name='us_rake', units='degrees', 
        dim=list(dim_unit_sources), missval=NA, 
        longname='rake_of_unit_source_in_degrees', 
        prec='float')
    us_slip_v = ncvar_def(name='us_slip', units='m', 
        dim=list(dim_unit_sources), missval=NA, 
        longname='slip_on_unit_source', 
        prec='float')
    us_length_v = ncvar_def(name='us_length', units='km', 
        dim=list(dim_unit_sources), missval=NA, 
        longname='length_of_unit_source', 
        prec='float')
    us_width_v = ncvar_def(name='us_width', units='km', 
        dim=list(dim_unit_sources), missval=NA, 
        longname='width_of_unit_source', 
        prec='float')
    us_downdip_number_v = ncvar_def(name='us_downdip_number', units='', 
        dim=list(dim_unit_sources), missval=NA, 
        longname='down_dip_index_of_unit_source', 
        prec='integer')
    us_alongstrike_number_v = ncvar_def(name='us_alongstrike_number', units='', 
        dim=list(dim_unit_sources), missval=NA, 
        longname='alongstrike_index_of_unit_source', 
        prec='integer')
    us_subfault_number_v = ncvar_def(name='us_subfault_number', units='', 
        dim=list(dim_unit_sources), missval=NA, 
        longname='subfault_index_of_unit_source', 
        prec='integer')
    us_max_depth_v = ncvar_def(name='us_max_depth', units='km', 
        dim=list(dim_unit_sources), missval=NA, 
        longname='max_depth_of_unit_source', prec='float')
    us_initial_condition_v = ncvar_def(name='us_initial_condition_file',
        units='', dim=list(dim_unit_sources), missval=NA, 
        longname='filename_containing_initial_surface_displacement_for_unit_source',
        prec='char')
    us_tide_gauge_file_v = ncvar_def(name='us_tide_gauge_file',
        units='', dim=list(dim_unit_sources), missval=NA, 
        longname='filename_containing_tide_gauge_time_series_for_unit_source',
        prec='char')

    # Store in netcdf variable list
    all_nc_var = c(all_nc_var, list(us_lon_c_v, us_lat_c_v, us_depth_v, 
        us_strike_v, us_dip_v, us_rake_v, us_slip_v, us_length_v, us_width_v,
        us_downdip_number_v, us_alongstrike_number_v, us_subfault_number_v,
        us_max_depth_v, us_initial_condition_v, us_tide_gauge_file_v) )


    #
    # Make file
    #
    sourcename_dot_nc = paste0(sourcename, '.nc')
    output_file_name = paste0(
        #Name of all_eq_events_file, with the 'sourcename.nc' cut off
        gsub(sourcename_dot_nc, '', all_eq_events_file, ),
        # add on _tsunami_sourcename.nc
        '_tsunami_', sourcename_dot_nc)

    output_nc_file = nc_create(output_file_name, vars=all_nc_var)

    #
    # Add global attributes
    #
    ncatt_put(output_nc_file, varid=0, attname='earthquake_events_file',
        attval=all_eq_events_file, prec='char')
    ncatt_put(output_nc_file, varid=0, attname='unit_source_statistics_file',
        attval=unit_source_statistics_file, prec='char')
    ncatt_put(output_nc_file, varid=0, attname='source_zone_name',
        attval=source_zone_name)

    #
    # Add variables: gauge locations
    #
    ncvar_put(output_nc_file, gauge_lon_v, gauge_locations$lon) 
    ncvar_put(output_nc_file, gauge_lat_v, gauge_locations$lat) 
    ncvar_put(output_nc_file, gauge_elev_v, gauge_locations$elev) 
    ncvar_put(output_nc_file, gauge_id_v, gauge_locations$gaugeID) 

    #
    # Add gauge summary statistics
    #
    ncvar_put(output_nc_file, gauge_event_max_stage_v, gauge_event_max_stage)
    ncvar_put(output_nc_file, gauge_event_reference_period_v, gauge_event_reference_period)
    ncvar_put(output_nc_file, gauge_event_peak_to_trough_v, gauge_event_peak_to_trough)
    ncvar_put(output_nc_file, gauge_event_arrival_time_v, gauge_event_arrival_time)
    ncvar_put(output_nc_file, gauge_event_initial_stage_v, gauge_event_initial_stage)
   
    #
    # Add event summary statistics 
    #
    ncvar_put(output_nc_file, event_area_v, all_eq_events$area)
    ncvar_put(output_nc_file, event_mean_length_v, all_eq_events$mean_length)
    ncvar_put(output_nc_file, event_mean_width_v, all_eq_events$mean_width)
    ncvar_put(output_nc_file, event_slip_v, all_eq_events$slip)
    ncvar_put(output_nc_file, event_Mw_v, all_eq_events$Mw)
    ncvar_put(output_nc_file, event_mean_depth_v, all_eq_events$mean_depth)
    ncvar_put(output_nc_file, event_max_depth_v, all_eq_events$max_depth)
    ncvar_put(output_nc_file, event_index_string_v, all_eq_events$event_index_string)
    ncvar_put(output_nc_file, event_sourcename_v, all_eq_events$sourcename)
    ncvar_put(output_nc_file, event_rate_annual_v, all_eq_events$rate_annual)
    ncvar_put(output_nc_file, event_rate_annual_upper_ci_v, all_eq_events$rate_annual_upper_ci)
    ncvar_put(output_nc_file, event_rate_annual_lower_ci_v, all_eq_events$rate_annual_lower_ci)

    #
    # Add unit source summary statistics 
    #
    ncvar_put(output_nc_file, us_lon_c_v, unit_source_statistics$lon_c)
    ncvar_put(output_nc_file, us_lat_c_v, unit_source_statistics$lat_c)
    ncvar_put(output_nc_file, us_depth_v, unit_source_statistics$depth)
    ncvar_put(output_nc_file, us_strike_v, unit_source_statistics$strike)
    ncvar_put(output_nc_file, us_dip_v, unit_source_statistics$dip)
    ncvar_put(output_nc_file, us_rake_v, unit_source_statistics$rake)
    ncvar_put(output_nc_file, us_slip_v, unit_source_statistics$slip)
    ncvar_put(output_nc_file, us_length_v, unit_source_statistics$length)
    ncvar_put(output_nc_file, us_width_v, unit_source_statistics$width)
    ncvar_put(output_nc_file, us_downdip_number_v, unit_source_statistics$downdip_number)
    ncvar_put(output_nc_file, us_alongstrike_number_v, unit_source_statistics$alongstrike_number)
    ncvar_put(output_nc_file, us_subfault_number_v, unit_source_statistics$subfault_number)
    ncvar_put(output_nc_file, us_max_depth_v, unit_source_statistics$max_depth)
    ncvar_put(output_nc_file, us_initial_condition_v, unit_source_statistics$initial_condition)
    ncvar_put(output_nc_file, us_tide_gauge_file_v, unit_source_statistics$tide_gauge_file)

    # Finish and flush to disk
    nc_close(output_nc_file)
    return(invisible(output_file_name))
}

