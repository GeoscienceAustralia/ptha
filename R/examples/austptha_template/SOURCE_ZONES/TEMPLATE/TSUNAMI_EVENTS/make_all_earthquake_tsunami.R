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

# Files to read
source_zone_name = basename(dirname(getwd()))
earthquake_events_file = paste0('all_uniform_slip_earthquake_events_', source_zone_name, '.nc')
unit_source_statistics_file = paste0('unit_source_statistics_', source_zone_name, '.nc')

# Run parameters
msl = 0.0 # Pre-tsunami initial condition for linear model. Gauges with elev>msl are inactive
lat_range = c(-72 + 2/60, 65 - 2/60) # Latitude range where gauges should be taken
lon_range = c(-40, 320) # Longitude range where gauges can be taken

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
full_unit_sources_mem = real_bytes * nvar * length(gauge_times) * ngauges * nunit_sources /(1024 * 1024)
station_chunk_size = floor( length(gauge_locations[,1]) * min(1, memory_for_unit_sources/(full_unit_sources_mem*mc_cores)) )
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

    output = matrix(NA, nrow=dim(flow_data)[1], ncol=4)

    for(i in 1:dim(flow_data)[1]){
        stages = flow_data[i,,1]
        rst = range(stages)
        # Peak stage
        output[i,1] = max(rst)
        # Reference period
        # NOTE: This period is approximate in a number of ways -- consider revising
        output[i,2] = rptha::zero_crossing_period(stages, dt=mean_dt)
        # Peak_to_trough
        output[i,3] = diff(rst)
        # Arrival time
        # Suppose > 1mm
        suppressWarnings({kk = min(which(abs(stages) > 1.0e-03))})
        if(is.finite(kk)){
            output[i,4] = times[kk]
        }

    }

    return(output)

}

## Loop over chunks of gauges, computing the statistics as we go
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

#
# Parallel version of the above loop
#
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



# Pack the statistics into the output arrays
for(i in 1:length(gauge_chunks_list)){
    gcl = gauge_chunks_list[[i]]

    # Vector with 1 at 'valid gauges', 0 elsewhere
    valid_gauges = (
        gauge_locations$lat[gcl] < (65-2/60) &
        gauge_locations$lat[gcl] > (-72+2/60) &
        gauge_locations$elev[gcl] < msl)

    # Set 'invalid gauges' values to zero
    for(j in 1:nevents){
        gauge_event_max_stage[j, gcl ] = modelled_flow_store[[i]][[j]][,1] * valid_gauges
        gauge_event_reference_period[j, gcl ] = modelled_flow_store[[i]][[j]][,2] * valid_gauges
        gauge_event_peak_to_trough[j, gcl ] = modelled_flow_store[[i]][[j]][,3] * valid_gauges
        gauge_event_arrival_time[j, gcl ] = modelled_flow_store[[i]][[j]][,4] * valid_gauges
    }
}


