# Record start time, so we can report how long the script takes
time_start = Sys.time()
print(paste0('Time is: ', time_start))

# Key libraries
suppressPackageStartupMessages(library(rptha))
suppressPackageStartupMessages(library(parallel))
source('sum_tsunami_unit_sources.R')
#
# Get local parameters code
#
config_env = new.env()
source('config.R', local=config_env)


###############################################################################

#
# Input parameters
#

command_arguments = commandArgs(trailingOnly=TRUE)

# Set memory we will allow for unit-sources tsunami, in MB. 
# It's a rough estimator -- should be a fraction of 1-nodes memory [e.g. using
# 10 GB on a 32GB machine seemed to work ok]
memory_for_unit_sources = config_env$memory_for_unit_sources  

# Number of cores in parallel [shared memory only]
mc_cores = config_env$mc_cores

# Run parameters

# Pre-tsunami initial condition for linear model. Gauges with elev>msl are
# inactive
msl = config_env$msl 

# Latitude range where gauges should be taken
lat_range = config_env$lat_range 

# Record the wave arrival time as the time at which abs(stage) exceeds this
# threshold
stage_threshold_for_arrival_time = config_env$stage_threshold_for_arrival_time 

# Files to read
source_zone_name = basename(dirname(getwd()))
unit_source_statistics_file = paste0('unit_source_statistics_', 
    source_zone_name, '.nc')
# Either uniform or stochastic or variable_uniform slip
# ... uniform slip is default
earthquake_events_file = paste0('all_uniform_slip_earthquake_events_', 
    source_zone_name, '.nc')
# .... but replace with stochastic slip if an argument -stochastic_slip was
# passed to R on startup, or request -variable_uniform_slip 
stochastic_slip = FALSE
if(length(command_arguments) > 0){

    if(any(grepl('-stochastic_slip', command_arguments))){
        earthquake_events_file = paste0('all_stochastic_slip_earthquake_events_', 
            source_zone_name, '.nc')
        stochastic_slip = TRUE
    }

    # We can also treat 'variable uniform slip' cases as stochastic slip
    # This is simple, though perhaps not the most efficient approach
    if(any(grepl('-variable_uniform_slip', command_arguments))){
        earthquake_events_file = paste0(
            'all_variable_uniform_slip_earthquake_events_', 
            source_zone_name, '.nc')
        stochastic_slip = TRUE
    }
}

# Setup output file name
sourcename_dot_nc = paste0(source_zone_name, '.nc')
output_file_name = paste0(
    #Name of earthquake_events_file, with the 'sourcename.nc' cut off
    gsub(sourcename_dot_nc, '', earthquake_events_file),
    # add on tsunami_sourcename.nc 
    'tsunami_', sourcename_dot_nc)



# Do we only make the file (if -make_file_only was passed), or do we run the
# computation (default)? 
make_file_only = FALSE
if(length(command_arguments) > 0){
    if(any(grepl('-make_file_only', command_arguments))){
        make_file_only=TRUE
        print('Making output file with empty flow variables')
        # Remove any RDS-for-chunk files that might be left over from previous
        # runs. Note this 'glob' corresponds to the 'image_file' created below
        # if save_as_RDS=TRUE
        old_RDS_files = Sys.glob(paste0(config_env$tmp_RDS_dir, '/',
            gsub('.nc', '', output_file_name), '_RDS_*.RDS'))
        unlink(old_RDS_files)

    }
}

# If -only_update_non_flow_variables was passed as a commandline argument,
# then we assume the output_file_name already exists, but that the non-flow
# variables require an update
only_update_non_flow_variables = FALSE
if(length(command_arguments) > 0){
    if(any(grepl('-only_update_non_flow_variables', command_arguments))){
        only_update_non_flow_variables = TRUE
        stopifnot(file.exists(output_file_name))
    }
}

# Check whether we only run a subset of events (this happens if commandline
# arguments of the form " --subset 4 10" where passed. Here '4' is this_subset,
# and 10 is number_of_subsets.) Subsetting gives us a straightforward way to run
# the code in parallel
subset_only=FALSE
if(any(grepl('-subset', command_arguments))){

    subset_only=TRUE

    command_ind = which(grepl('-subset', command_arguments))
    stopifnot(length(command_ind)== 1)

    # Subset will be followed by two integers, giving the 
    # index of 'this_subset', and the total number of subsets
    this_subset = as.numeric(command_arguments[command_ind+1])
    number_of_subsets = as.numeric(command_arguments[command_ind+2])
    stopifnot(this_subset <= number_of_subsets)

    print(paste0('running subset ', this_subset, '/', number_of_subsets))

    # Better check that the output file exists
    if(!file.exists(output_file_name)){
        stop(paste0('Could not find output file ', output_file_name))
    }

    # Check for save_as_RDS option
    save_as_RDS = FALSE
    if(any(grepl('-save_as_RDS', command_arguments))){
        save_as_RDS = TRUE
    }

}

#
# End input
#
###############################################################################


# Read unit source statistics & uniform slip earthquake events
unit_source_statistics = read_table_from_netcdf(unit_source_statistics_file)
all_eq_events = read_table_from_netcdf(earthquake_events_file)
# Find max number of characters required to store columns -- since this should
# be passed as a lower bound to the output netcdf file
max_nchar_unit_source_statistics = max(unlist(lapply(unit_source_statistics, 
    f<-function(x) max(nchar(x)*is.character(x)))))
max_nchar_all_eq_events = max(unlist(lapply(all_eq_events, 
    f<-function(x) max(nchar(x)*is.character(x)))))


# Potentially only look at a subset of events
my_events = 1:length(all_eq_events[,1])
if(subset_only){
    # Update my_events and all_eq_events
    my_events = splitIndices(length(all_eq_events[,1]), 
        number_of_subsets)[[this_subset]]
    all_eq_events = all_eq_events[my_events,]
    gc()
}
# Ensure a contiguous chunk was selected
stopifnot(all( range(diff(my_events)) == 1))

# Read gauges
gauge_locations = get_netcdf_gauge_locations(
    unit_source_statistics$tide_gauge_file[1])
# Read times at which gauge values are stored -- always constant thanks to the
# model setup
gauge_times = get_netcdf_gauge_output_times(
    unit_source_statistics$tide_gauge_file[1])

#
# Partition the work among processes [shared memory parallel]
#
nevents = length(all_eq_events[,1])
ngauges = length(gauge_locations[,1])
nunit_sources = length(unit_source_statistics[,1])
real_bytes = 8  # size of an R real number
nvar = 3 # [stage, uh, vh]
full_unit_sources_mem = real_bytes * nvar * length(gauge_times) * ngauges * 
    nunit_sources /(1024 * 1024)
# Make the station_chunk_size so that when spread over all cores, it 
# takes up memory = memory_for_unit_sources.
# My 'model' of how this code uses memory is very approximate and
# under-estimates, so you should set 'memory_for_unit_sources' to be say 1/3 of
# the available memory.
station_chunk_size = floor( ngauges * 
    min(1, memory_for_unit_sources/(full_unit_sources_mem*mc_cores)) )
gauge_chunks_list = splitIndices(ngauges, ceiling(ngauges/station_chunk_size))

#
# Given an array with [gauge, time-slice, [stg] ], extract flow summary
# statistics for each gauge. Pass this to the function that makes the tsunami
# events, so we don't have to store full time-series [too much memory for
# many gauges and events]
#
local_summary_function<-function(flow_data){

    output = gauge_statistics_simple(gauge_times, flow_data, 
        stage_threshold_for_arrival_time)
    return(output)

}

# Function to compute tsunami events for a 'chunk' of gauges. Used to farm
# out the work in parallel
parfun<-function(gcl){
    make_tsunami_event_from_unit_sources(
        earthquake_events = all_eq_events,
        unit_source_statistics = unit_source_statistics,
        unit_source_flow_files = unit_source_statistics$tide_gauge_file,
        indices_of_subset = gcl,
        summary_function = local_summary_function,
        msl = msl,
        all_flow_variables=FALSE)
}

if(!make_file_only){
    # Here we avoid use of mclapply, which seems to leave un-stopped worker
    # nodes on NCI
    cl = makeForkCluster(nnodes=mc_cores)
    modelled_flow_store = parLapply(cl=cl, X=gauge_chunks_list, fun=parfun)
    stopCluster(cl)
}

#
# Pack outputs
#

# Make a 'missing_data' flag which is definitely interpreted as real by netcdf
nul_r = config_env$null_double 

if(!only_update_non_flow_variables){
    # Matrix for max stage, with as many rows as events, and as many columns as
    # gauges.
    # Even if (make_file_only==TRUE), we use this to create the netcdf output
    gauge_event_max_stage = matrix(nul_r, nrow=nevents, ncol=ngauges)

}else{
    gauge_event_max_stage = NULL
    
}

if(make_file_only | only_update_non_flow_variables){

    # Dummy values to save memory
    gauge_event_reference_period = NULL
    gauge_event_peak_to_trough = NULL
    gauge_event_arrival_time = NULL
    gauge_event_initial_stage = NULL

}else{

    # Period
    gauge_event_reference_period = matrix(nul_r, nrow=nevents, ncol=ngauges)
    # Peak-to-trough
    gauge_event_peak_to_trough = matrix(nul_r, nrow=nevents, ncol=ngauges)
    # Arrival time
    gauge_event_arrival_time = matrix(nul_r, nrow=nevents, ncol=ngauges)
    # Initial stage
    gauge_event_initial_stage = matrix(nul_r, nrow=nevents, ncol=ngauges)

    # Pack the statistics into the output arrays
    for(i in 1:length(gauge_chunks_list)){
        gcl = gauge_chunks_list[[i]]

        # Vector with 1 at 'valid gauges', 0 elsewhere
        valid_gauges = (
            gauge_locations$lat[gcl] < lat_range[2] &
            gauge_locations$lat[gcl] > lat_range[1] &
            gauge_locations$elev[gcl] < msl)

        valid_gauges_m = (c(NA, 1)[valid_gauges+1])

        # Set 'invalid gauges' values to NA
        for(j in 1:nevents){
            gauge_event_max_stage[j, gcl ] = 
                modelled_flow_store[[i]][[j]][,1] * valid_gauges_m 
            gauge_event_reference_period[j, gcl ] = 
                modelled_flow_store[[i]][[j]][,2] * valid_gauges_m
            gauge_event_peak_to_trough[j, gcl ] = 
                modelled_flow_store[[i]][[j]][,3] * valid_gauges_m
            gauge_event_arrival_time[j, gcl ] = 
                modelled_flow_store[[i]][[j]][,4] * valid_gauges_m
            gauge_event_initial_stage[j, gcl ] = 
                modelled_flow_store[[i]][[j]][,5] * valid_gauges_m
        }
    }

    # Forcibly free some memory
    rm(modelled_flow_store, valid_gauges, valid_gauges_m)
    gc()
}
###############################################################################
#
# Write all the relevant information to netcdf
#
###############################################################################

write_all_source_zone_tsunami_statistics_to_netcdf<-function(
    source_zone_name,
    all_eq_events,
    earthquake_events_file,
    unit_source_statistics,
    unit_source_statistics_file,
    gauge_locations,
    gauge_event_max_stage,
    gauge_event_reference_period,
    gauge_event_peak_to_trough,
    gauge_event_arrival_time,
    gauge_event_initial_stage,
    stage_threshold_for_arrival_time,
    max_nchar_all_eq_events,
    max_nchar_unit_source_statistics,
    output_file_name,
    only_update_non_flow_variables){

    library(ncdf4)

    if(!only_update_non_flow_variables){

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
        charlen_initial_cond_string = max(nchar(
            unit_source_statistics$initial_condition_file))
        charlen_tide_gauge_string = max(nchar(
            unit_source_statistics$tide_gauge_file))

        dim_char_size = max(c(charlen_sourcename, charlen_event_index_string, 
            charlen_initial_cond_string, charlen_tide_gauge_string,
            max_nchar_all_eq_events, max_nchar_unit_source_statistics))

        dim_nchar = ncdim_def(name='max_nchar', units='', vals=1:dim_char_size, 
            unlim=FALSE,
            longname='integer index corresponding to the maximum number of string characters')


        #
        # Create netcdf variable for gauges
        #
        gauge_lon_v = ncvar_def(name='lon', units='degrees_east', 
            dim=list(dim_station), missval=NA, longname='station_longitude', 
            prec='float')
        gauge_lat_v = ncvar_def(name='lat', units='degrees_north', 
            dim=list(dim_station), missval=NA, longname='station_latitude', 
            prec='float')
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
        # These statistics vary depending on whether we do stochastic_slip or
        # uniform_slip
        #
        if(stochastic_slip == FALSE){
            #
            # Uniform slip case
            #
            event_area_v = ncvar_def(name='event_area', units='km^2', 
                dim=list(dim_event), missval=NA, 
                longname='rupture_area_of_earthquake_event', prec='float')

            event_mean_length_v = ncvar_def(name='event_mean_length', units='km', 
                dim=list(dim_event), missval=NA, 
                longname='mean_rupture_length_of_earthquake_event', prec='float')

            event_mean_width_v = ncvar_def(name='event_mean_width', units='km', 
                dim=list(dim_event), missval=NA, 
                longname='mean_rupture_width_of_earthquake_event', prec='float')

            event_slip_v = ncvar_def(name='event_slip', units='m', 
                dim=list(dim_event), missval=NA, 
                longname='slip_of_earthquake_event', prec='float')

            event_Mw_v = ncvar_def(name='event_Mw', units='', dim=list(dim_event),
                missval=NA, longname='moment_magnitude_of_earthquake_event', 
                prec='float')

            event_mean_depth_v = ncvar_def(name='event_mean_depth', units='km', 
                dim=list(dim_event), missval=NA, 
                longname='mean_rupture_depth_of_earthquake_event', prec='float') 

            event_max_depth_v = ncvar_def(name='event_max_depth', units='km', 
                dim=list(dim_event), missval=NA, 
                longname='max_rupture_depth_of_earthquake_event', prec='float') 

            event_index_string_v = ncvar_def(name='event_index_string', units='', 
                dim=list(dim_nchar, dim_event), missval=NULL, 
                longname='indices_of_unit_sources_included_in_earthquake_event', 
                prec='char')
            event_sourcename_v = ncvar_def(name='event_sourcename', units='', 
                dim=list(dim_nchar, dim_event), missval=NULL, 
                longname='source_zone_name_of_earthquake_event', prec='char')

            event_rate_annual_v = ncvar_def(name='event_rate_annual', 
                units='events per year', dim=list(dim_event),
                missval=NA, longname='mean_rate_of_earthquake_event', prec='float') 
            event_rate_annual_upper_ci_v = ncvar_def(
                name='event_rate_annual_upper_ci', 
                units='events per year', dim=list(dim_event),
                missval=NA, 
                longname='upper_credible_interval_for_rate_of_earthquake_event', 
                prec='float') 
            event_rate_annual_lower_ci_v = ncvar_def(
                name='event_rate_annual_lower_ci', 
                units='events per year', dim=list(dim_event),
                missval=NA, 
                longname='lower_credible_interval_for_rate_of_earthquake_event', 
                prec='float')

            # Append to list of variables, for making netcdf
            all_nc_var = c(all_nc_var, list(event_area_v, event_mean_length_v, 
                event_mean_width_v, event_slip_v, event_Mw_v, event_mean_depth_v, 
                event_max_depth_v, event_index_string_v, event_sourcename_v, 
                event_rate_annual_v, event_rate_annual_upper_ci_v, 
                event_rate_annual_lower_ci_v) )
        }else{
            #
            # Stochastic slip case
            #
            event_Mw_v = ncvar_def(name='event_Mw', units='', dim=list(dim_event),
                missval=NA, longname='moment_magnitude_of_earthquake_event', 
                prec='float')
            event_target_lon_v = ncvar_def(name='event_target_lon', 
                units='degrees_east', dim=list(dim_event),
                missval=NA, longname='longitude_near_earthquake_event', 
                prec='float')
            event_target_lat_v = ncvar_def(name='event_target_lat', 
                units='degrees_north', dim=list(dim_event),
                missval=NA, longname='latitude_near_earthquake_event', prec='float')

            event_peak_slip_downdip_ind_v = ncvar_def(
                name='event_peak_slip_downdip_index', units='', 
                dim=list(dim_event), missval=NULL, 
                longname='down-dip_index_of_the_unit_source_with_peak_slip',
                prec='integer')
            event_peak_slip_alongstrike_ind_v = ncvar_def(
                name='event_peak_slip_alongstrike_index', units='', 
                dim=list(dim_event), missval=NULL, 
                longname='along-strike_index_of_the_unit_source_with_peak_slip',
                prec='integer')
            event_uniform_event_row_v = ncvar_def(
                name='event_uniform_event_row', units='', 
                dim=list(dim_event), missval=NULL, 
                longname='row_index_of_the_uniform_slip_event_used_to_define_location_and_magnitude',
                prec='integer')
            
            event_sourcename_v = ncvar_def(name='event_sourcename', units='', 
                dim=list(dim_nchar, dim_event),
                missval=NULL, longname='source_zone_name_of_earthquake_event', 
                prec='char')

            event_index_string_v = ncvar_def(name='event_index_string', units='', 
                dim=list(dim_nchar, dim_event),
                missval=NULL, 
                longname='indices_of_unit_sources_included_in_earthquake_event', 
                prec='char')
            event_slip_string_v = ncvar_def(name='event_slip_string', units='', 
                dim=list(dim_nchar, dim_event), missval=NULL, 
                longname='slip_(m)_of_unit_sources_included_in_earthquake_event_with_underscore_separator', 
                prec='char')

            event_rate_annual_v = ncvar_def(name='event_rate_annual', 
                units='events per year', dim=list(dim_event),
                missval=NA, longname='mean_rate_of_earthquake_event', prec='float') 
            event_rate_annual_upper_ci_v = ncvar_def(
                name='event_rate_annual_upper_ci', 
                units='events per year', dim=list(dim_event),
                missval=NA, 
                longname='upper_credible_interval_for_rate_of_earthquake_event', 
                prec='float') 
            event_rate_annual_lower_ci_v = ncvar_def(
                name='event_rate_annual_lower_ci', 
                units='events per year', dim=list(dim_event),
                missval=NA, 
                longname='lower_credible_interval_for_rate_of_earthquake_event', 
                prec='float')

            all_nc_var = c(all_nc_var, list(event_Mw_v, event_target_lon_v, 
                event_target_lat_v, event_peak_slip_downdip_ind_v, 
                event_peak_slip_alongstrike_ind_v, event_uniform_event_row_v,
                event_sourcename_v, event_index_string_v, event_slip_string_v,
                event_rate_annual_v, event_rate_annual_upper_ci_v,
                event_rate_annual_lower_ci_v))

        }

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
            dim=list(dim_unit_sources), missval=NULL, 
            longname='down_dip_index_of_unit_source', 
            prec='integer')
        us_alongstrike_number_v = ncvar_def(name='us_alongstrike_number', units='', 
            dim=list(dim_unit_sources), missval=NULL, 
            longname='alongstrike_index_of_unit_source', 
            prec='integer')
        us_subfault_number_v = ncvar_def(name='us_subfault_number', units='', 
            dim=list(dim_unit_sources), missval=NULL, 
            longname='subfault_index_of_unit_source', 
            prec='integer')
        us_max_depth_v = ncvar_def(name='us_max_depth', units='km', 
            dim=list(dim_unit_sources), missval=NA, 
            longname='max_depth_of_unit_source', prec='float')
        us_initial_condition_v = ncvar_def(name='us_initial_condition_file',
            units='', dim=list(dim_nchar, dim_unit_sources), missval=NULL, 
            longname='filename_containing_initial_surface_displacement_for_unit_source',
            prec='char')
        us_tide_gauge_file_v = ncvar_def(name='us_tide_gauge_file',
            units='', dim=list(dim_nchar, dim_unit_sources), missval=NULL, 
            longname='filename_containing_tide_gauge_time_series_for_unit_source',
            prec='char')

        # Store in netcdf variable list
        all_nc_var = c(all_nc_var, list(us_lon_c_v, us_lat_c_v, us_depth_v, 
            us_strike_v, us_dip_v, us_rake_v, us_slip_v, us_length_v, us_width_v,
            us_downdip_number_v, us_alongstrike_number_v, us_subfault_number_v,
            us_max_depth_v, us_initial_condition_v, us_tide_gauge_file_v) )


        output_nc_file = nc_create(output_file_name, vars=all_nc_var)

    }else{
        # Assume the file already exists, and just update some variables
        output_nc_file = nc_open(output_file_name)
    }

    #
    # Add global attributes
    #
    ncatt_put(output_nc_file, varid=0, attname='earthquake_events_file',
        attval=normalizePath(earthquake_events_file), prec='text')
    ncatt_put(output_nc_file, varid=0, attname='unit_source_statistics_file',
        attval=normalizePath(unit_source_statistics_file), prec='text')
    ncatt_put(output_nc_file, varid=0, attname='source_zone_name',
        attval=source_zone_name, prec='text')
    ncatt_put(output_nc_file, varid=0, attname='parent_script_name',
        attval=parent_script_name(), prec='text')

    #
    # Add variables: gauge locations
    #
    ncvar_put(output_nc_file, 'lon', gauge_locations$lon) ; gc()
    ncvar_put(output_nc_file, 'lat', gauge_locations$lat) ; gc()
    ncvar_put(output_nc_file, 'elev', gauge_locations$elev) ; gc()
    ncvar_put(output_nc_file, 'gaugeID', gauge_locations$gaugeID) ; gc()

    if(!only_update_non_flow_variables){
        #
        # Add gauge summary statistics
        #
        if(make_file_only){
            # Use gauge_event_max_stage to represent all output statistics
            ncvar_put(output_nc_file, 'max_stage', gauge_event_max_stage)
            gc()
            ncvar_put(output_nc_file, 'period', gauge_event_max_stage)
            gc()
            ncvar_put(output_nc_file, 'stage_range', gauge_event_max_stage)
            gc()
            ncvar_put(output_nc_file, 'arrival_time', gauge_event_max_stage)
            gc()
            ncvar_put(output_nc_file, 'initial_stage', gauge_event_max_stage)
            gc()

        }else{
            ncvar_put(output_nc_file, 'max_stage', gauge_event_max_stage)
            gc()
            ncvar_put(output_nc_file, 'period', gauge_event_reference_period)
            gc()
            ncvar_put(output_nc_file, 'stage_range', gauge_event_peak_to_trough)
            gc()
            ncvar_put(output_nc_file, 'arrival_time', gauge_event_arrival_time)
            gc()
            ncvar_put(output_nc_file, 'initial_stage', gauge_event_initial_stage)
            gc()
        }
    }
   
    #
    # Add event summary statistics 
    #
    if(stochastic_slip == FALSE){
        # Uniform slip
        ncvar_put(output_nc_file, 'event_area', all_eq_events$area); gc()
        ncvar_put(output_nc_file, 'event_mean_length', all_eq_events$mean_length); gc()
        ncvar_put(output_nc_file, 'event_mean_width', all_eq_events$mean_width); gc()
        ncvar_put(output_nc_file, 'event_slip', all_eq_events$slip); gc()
        ncvar_put(output_nc_file, 'event_Mw', all_eq_events$Mw); gc()
        ncvar_put(output_nc_file, 'event_mean_depth', all_eq_events$mean_depth); gc()
        ncvar_put(output_nc_file, 'event_max_depth', all_eq_events$max_depth); gc()
        ncvar_put(output_nc_file, 'event_index_string', all_eq_events$event_index_string); gc()
        ncvar_put(output_nc_file, 'event_sourcename', all_eq_events$sourcename); gc()
        ncvar_put(output_nc_file, 'event_rate_annual', all_eq_events$rate_annual); gc()
        ncvar_put(output_nc_file, 'event_rate_annual_upper_ci', all_eq_events$rate_annual_upper_ci); gc()
        ncvar_put(output_nc_file, 'event_rate_annual_lower_ci', all_eq_events$rate_annual_lower_ci); gc()
    }else{
        # Stochastic slip
        ncvar_put(output_nc_file, 'event_Mw', all_eq_events$Mw); gc()
        ncvar_put(output_nc_file, 'event_target_lon', all_eq_events$target_lon); gc()
        ncvar_put(output_nc_file, 'event_target_lat', all_eq_events$target_lat); gc()
        ncvar_put(output_nc_file, 'event_peak_slip_downdip_index', all_eq_events$peak_slip_downdip_ind); gc()
        ncvar_put(output_nc_file, 'event_peak_slip_alongstrike_index', all_eq_events$peak_slip_alongstrike_ind); gc()
        ncvar_put(output_nc_file, 'event_sourcename', all_eq_events$sourcename); gc()
        ncvar_put(output_nc_file, 'event_uniform_event_row', all_eq_events$uniform_event_row); gc()
        ncvar_put(output_nc_file, 'event_rate_annual', all_eq_events$rate_annual); gc()
        ncvar_put(output_nc_file, 'event_rate_annual_upper_ci', all_eq_events$rate_annual_upper_ci); gc()
        ncvar_put(output_nc_file, 'event_rate_annual_lower_ci', all_eq_events$rate_annual_lower_ci); gc()
        ncvar_put(output_nc_file, 'event_index_string', all_eq_events$event_index_string); gc()
        ncvar_put(output_nc_file, 'event_slip_string', all_eq_events$event_slip_string); gc()
    }

    #
    # Add unit source summary statistics 
    #
    ncvar_put(output_nc_file, 'us_lon_c', unit_source_statistics$lon_c); gc()
    ncvar_put(output_nc_file, 'us_lat_c', unit_source_statistics$lat_c); gc()
    ncvar_put(output_nc_file, 'us_depth', unit_source_statistics$depth); gc()
    ncvar_put(output_nc_file, 'us_strike', unit_source_statistics$strike); gc()
    ncvar_put(output_nc_file, 'us_dip', unit_source_statistics$dip); gc()
    ncvar_put(output_nc_file, 'us_rake', unit_source_statistics$rake); gc()
    ncvar_put(output_nc_file, 'us_slip', unit_source_statistics$slip); gc()
    ncvar_put(output_nc_file, 'us_length', unit_source_statistics$length); gc()
    ncvar_put(output_nc_file, 'us_width', unit_source_statistics$width); gc()
    ncvar_put(output_nc_file, 'us_downdip_number', unit_source_statistics$downdip_number); gc()
    ncvar_put(output_nc_file, 'us_alongstrike_number', unit_source_statistics$alongstrike_number); gc()
    ncvar_put(output_nc_file, 'us_subfault_number', unit_source_statistics$subfault_number); gc()
    ncvar_put(output_nc_file, 'us_max_depth', unit_source_statistics$max_depth); gc()
    ncvar_put(output_nc_file, 'us_initial_condition_file', unit_source_statistics$initial_condition); gc()
    ncvar_put(output_nc_file, 'us_tide_gauge_file', unit_source_statistics$tide_gauge_file); gc()

    # Finish and flush to disk
    nc_close(output_nc_file)
    return(invisible(output_file_name))
}


if(make_file_only | (subset_only==FALSE) | only_update_non_flow_variables){
    #
    # Write it out
    #
    write_all_source_zone_tsunami_statistics_to_netcdf(
        source_zone_name,
        all_eq_events,
        earthquake_events_file,
        unit_source_statistics,
        unit_source_statistics_file,
        gauge_locations,
        gauge_event_max_stage,
        gauge_event_reference_period,
        gauge_event_peak_to_trough,
        gauge_event_arrival_time,
        gauge_event_initial_stage,
        stage_threshold_for_arrival_time,
        max_nchar_all_eq_events,
        max_nchar_unit_source_statistics,
        output_file_name,
        only_update_non_flow_variables)

}else{
    # Save the variables -- either direct to the netcdf, or to an individual RDS file

    if(save_as_RDS){
        #
        # Save R image, so we can write to netcdf later
        #
        # This is useful to avoid accidently having multiple images writing to
        # netcdf at once
        #
        print('Writing to RDS')

        # Make directory for R images
        tmp_RDS_dir = config_env$tmp_RDS_dir # 'R_images_tmp'
        dir.create(tmp_RDS_dir, showWarnings=FALSE)

        image_file = paste0(tmp_RDS_dir, '/', 
            gsub('.nc', '', output_file_name), '_RDS_', 
            this_subset, '_', number_of_subsets, '.RDS')
        save.image(image_file)

    }else{
        #
        # Update variables in the already-created netcdf_file
        #

        # Open the file for editing -- DO NOT DO THIS WITH MULTIPLE PROGRAMS AT ONCE
        output_nc_file = nc_open(output_file_name, readunlim=FALSE, write=TRUE)

        # Put each variable, only in the contiguous part of my_events
        ncvar_put(output_nc_file, output_nc_file$var$max_stage, gauge_event_max_stage, 
                  start=c(my_events[1], 1), count=c(length(my_events), -1))
        gc()

        ncvar_put(output_nc_file, output_nc_file$var$period, gauge_event_reference_period,
                  start=c(my_events[1], 1), count=c(length(my_events), -1))
        gc()

        ncvar_put(output_nc_file, output_nc_file$var$stage_range, gauge_event_peak_to_trough,
                  start=c(my_events[1], 1), count=c(length(my_events), -1))
        gc()

        ncvar_put(output_nc_file, output_nc_file$var$arrival_time, gauge_event_arrival_time,
                  start=c(my_events[1], 1), count=c(length(my_events), -1))
        gc()

        ncvar_put(output_nc_file, output_nc_file$var$initial_stage, gauge_event_initial_stage,
                  start=c(my_events[1], 1), count=c(length(my_events), -1))
        gc()

        nc_close(output_nc_file)
    }
}

# Some useful finishing information
time_end = Sys.time()
print(paste0('Ending time: ', time_end))
print(paste0('    nevents: ', length(my_events)))
print(paste0('That_took: ', format(time_end - time_start, units='s')))
print(gc())
