#
# Compute the minimum and scenario_average arrival times.
#
# In our SWALS calculations, the arrival time was defined as the first time at which
#   max-stage > 0.01 + background-sea-level
#
# The intention of the 'minimum' arrival times is to give an idea of how much
# time they should have before a tsunami arrives, so they can be out of the
# danger zone when the wave arrives.
#
# The intention of the scenario_average arrival times is to give an indication
# of how much later than the 'minimum arrival time' a tsunami might arrive. It
# will help emergency services think about how much additional time they might
# have, and what could plausibly be achieved.
# - We don't want to weight the average by the regular scenario rate, because then
#   few scenarios will contribute much (i.e. those with high rates and small
#   tsunamis).
# - A simple alternative is to average over all scenarios without regard to
#   their rate. This result will be sensitive to how we setup the importance
#   sampling. It will also be sensitive to how we define the arrival time.
#   Nonetheless it will give an indication
#

#
# INPUTS
#
domain_indices_to_process = seq(1, 488)
raster_output_dir = './nsw_mis_arrival_time_min_and_scenario_average'

#
# END INPUTS
#

suppressPackageStartupMessages(library(utils))
suppressPackageStartupMessages(library(rptha))
suppressPackageStartupMessages(library(terra))

# Probabilistic_inundation metadata: To link scenario rates to scenarios
asfm = new.env()
source('application_specific_file_metadata.R', local=asfm, chdir=TRUE)

# Get source zone metadata
source_zones = names(asfm$source_zone_modelled_tsunami_scenario_basedirs)
source_zones_metadata = lapply(asfm$source_zone_modelled_tsunami_scenario_basedirs, 
    asfm$get_scenario_metadata_from_md_base_dir)


process_domain_index<-function(domain_index, raster_tar_files, output_dir){
    
    # Useful to have one of the rasters ready as a template [to help with data
    # The raster depth file associated with each md_dir. There should only be one
    # per md_dir. All rasters must have the same extent and resolution.
    raster_files_one_domain = paste0('/vsitar/', raster_tar_files, 
        "/arrival_time_domain_", domain_index, ".tif")

    # Make space to store the output
    raster_template = rast(raster_files_one_domain[1])
    # Start with zeros everywhere (even NA cells)
    raster_template[is.na(raster_template)] = 0.0
    raster_template = setValues(raster_template, 0.0)

    # The calculations require storing the minimum, and two values for the average (numerator and denominator)
    HUGE_VALUE = 9e+50
    min_arrival_time = as.matrix(raster_template, wide=TRUE) + HUGE_VALUE
    sum_arrival_time = as.matrix(raster_template, wide=TRUE)
    arrival_time_nonNA_count = as.matrix(raster_template[[1]], wide=TRUE)

    for(i in 1:length(raster_files_one_domain)){
        # Convert raster_i to an array
        raster_i = rast(raster_files_one_domain[i])
        result_i = as.matrix(raster_i, wide=TRUE)

        stopifnot(all(dim(result_i) == dim(min_arrival_time)))

        is_na_result_i = is.na(result_i)
       
        # Count number of times each cell is not missing data (denominator for average calculation)
        arrival_time_nonNA_count = arrival_time_nonNA_count + (1 - is_na_result_i)
        # Numerator for average calculation
        result_i[is.na(result_i)] = 0
        sum_arrival_time = sum_arrival_time + result_i

        # Minimum arrival time
        result_i[result_i <= 0] = HUGE_VALUE
        min_arrival_time = pmin(min_arrival_time, result_i)

    }

    # Make minimum arrival time raster
    min_arrival_time_rast = raster_template
    min_arrival_time_rast = setValues(min_arrival_time_rast, min_arrival_time) 
    # Treat cells that didn't have an arrival time
    min_arrival_time_rast[min_arrival_time_rast > 0.99*HUGE_VALUE] = NA

    # Make the mean arrival time raster
    mean_arrival_time_rast = raster_template
    # Convert to average (avoiding division by zero)
    sum_arrival_time = sum_arrival_time/pmax(arrival_time_nonNA_count, 1e-06)
    # Treat cells that didn't have an arrival time
    sum_arrival_time[arrival_time_nonNA_count == 0] = NA  
    mean_arrival_time_rast = setValues(mean_arrival_time_rast, sum_arrival_time)

    #output = list(min_arrival_time_rast = min_arrival_time_rast,
    #              mean_arrival_time_rast = mean_arrival_time_rast,
    #              domain_index = domain_index)

    output_file = paste0(output_dir, '/minimum_arrival_time_domain_', domain_index, '.tif')
    writeRaster(min_arrival_time_rast, file=output_file, gdal=c('COMPRESS=DEFLATE'), overwrite=TRUE)
    output_file = paste0(output_dir, '/mean_arrival_time_domain_', domain_index, '.tif')
    writeRaster(mean_arrival_time_rast, file=output_file, gdal=c('COMPRESS=DEFLATE'), overwrite=TRUE)

    rm(raster_template, min_arrival_time, sum_arrival_time, 
        arrival_time_nonNA_count, raster_i, result_i, is_na_result_i, 
        min_arrival_time_rast, mean_arrival_time_rast)
    gc()
    return(0)
}

# Prevent failures in parallel from causing the entire computation to fail
try_process_domain_index<-function(domain_index, raster_tar_files, output_dir){
    try(process_domain_index(domain_index, raster_tar_files, output_dir))
}

#
# Run the calculations
#    

# Get the rasters and their rate information (albeit rate information not currently used)
#source_zones_rates_and_rasters = get_logic_tree_mean_rates_and_md_files_all_source_zones(
#    source_zones_scenarios, source_zones_metadata)

# Make output folders (one per source zone)
dir.create(raster_output_dir, showWarnings=FALSE, recursive=TRUE)
output_dir_sz = paste0(raster_output_dir, '/', names(source_zones_metadata))
for(i in 1:length(output_dir_sz)){
    dir.create(output_dir_sz[i], showWarnings=FALSE)
}

library(parallel)
MC_CORES = asfm$DEFAULT_MC_CORES
# Setup a cluster
local_cluster = makeCluster(MC_CORES)

# Load packages on cluster
ignore = clusterCall(local_cluster, fun=function(){ 
    suppressPackageStartupMessages(library(utils)); 
    suppressPackageStartupMessages(library(rptha)); 
    suppressPackageStartupMessages(library(terra)) })

# Copy over all variables
vars_to_export = ls(all=TRUE)
clusterExport(local_cluster, varlist=vars_to_export)

# Run the calculation for all source-zones
for(i in 1:length(source_zones_metadata)){
    sz = names(source_zones_metadata)[i]
    output_dir_i = output_dir_sz[i] 

    # Main parallel calculation
    all_raster_calculations = parLapplyLB(cl = local_cluster, 
        X = domain_indices_to_process, 
        fun = try_process_domain_index, 
        raster_tar_files = source_zones_metadata[[sz]]$raster_tar_files,
        output_dir = output_dir_i,
        chunk.size=1)

    # Check for errors, and report if needed
    any_failed = unlist(lapply(all_raster_calculations, function(x) is(x, 'try-error')))
    if(any(any_failed)){
        print(paste0('WARNING: calculations failed for the following domain indices on source_zone: ', sz))
        k = which(any_failed)
        print(domain_indices_to_process[k])
    }
}

# Shut down the cluster
stopCluster(local_cluster)
