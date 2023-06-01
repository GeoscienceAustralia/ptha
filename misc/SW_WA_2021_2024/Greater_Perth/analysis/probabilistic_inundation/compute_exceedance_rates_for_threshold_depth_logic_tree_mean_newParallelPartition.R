#
# Compute rasters with rate (depth > threshold).
# Can also do similar computations for max-stage, arrival-time, etc.
#
#library(terra)
library(utils)
suppressPackageStartupMessages(library(rptha))
suppressPackageStartupMessages(library(parallel))
source('exceedance_rate_raster_calculations.R')


usage_msg = paste0(
    'Usage is like this (e.g. to compute results for domains 1:4 using runs ', 
    'in ../../swals/OUTPUTS/ptha18_tonga_MSL0):\n',
    '    Rscript compute_exceedance_rates_for_threshold_depth_logic_tree_mean_newParallelPartition.R ptha18_tonga_MSL0 4')

#
# INPUTS
#

# To facilitate running multiple cases, it is easy to pass some input arguments
# via the commandline
input_args = commandArgs(trailingOnly=TRUE)
if(!(length(input_args) == 2 & is.finite(as.numeric(input_args[2])))) stop(usage_msg)

# The name of the folder [beneath swals/OUTPUTS/] where the runs for this
# series are located.
run_series_name = input_args[1] 
# The maximum domain index (corresponding to an integer in the tif filename).
# We do calculations for domain 1, 2, 3, ... max_domain_index
max_domain_index = as.numeric(input_args[2]) 

# Number of cores -- the exceedance-rate computation is distributed over these.
MC_CORES = 48 
if(MC_CORES > detectCores()){
    stop(' MC_CORES exceeds the number of cores on your machine. Bad idea.')
}

# Rasters that we operate on must have filename starting with this.
# By changing this variable we can adapt calculations to:
#     max_stage ("max_stage_domain_")
#     arrival_time ("arrival_time_domain_")
# among other things.
raster_name_start = "depth_as_max_stage_minus_elevation0_domain_"

# Source-zone specific info
source_info = basename(run_series_name)
if(source_info == 'random_outerrisesunda'){
    # Key parameters for outerrisesunda source

    source_zone = 'outerrisesunda'
    scenario_base = '../../sources/hazard/random_outerrisesunda/'

    # Files with scenario row indices and other metadata
    scenario_data = paste0(scenario_base, 
        c('random_scenarios_outerrisesunda_unsegmented_HS.csv')
        )

    # Names for each source representation in scenario_data [unsegmented, and
    # various segments]
    names_scenario_data = gsub('.csv', '',
        gsub('random_scenarios_outerrisesunda_', '', basename(scenario_data),
             fixed=TRUE),
        fixed=TRUE)

    # Tarred multidomain directories for all SWALS model runs, for every scenario.
    # These do NOT need to be ordered in the same was as the scenario_data.
    tarred_multidomain_dirs = Sys.glob(paste0('../../swals/OUTPUTS/', run_series_name, 
        '/ptha18_random_scenarios_outerrisesunda_row_*/RUN*'))

}else if(source_info == 'random_sunda2'){
    # Key parameters for sunda2 source

    source_zone = 'sunda2'
    scenario_base = '../../sources/hazard/random_sunda2/'

    # Files with scenario row indices and other metadata
    scenario_data = paste0(scenario_base, 
        c('random_scenarios_sunda2_andaman_segment_HS.csv',
          'random_scenarios_sunda2_arakan_segment_HS.csv',
          'random_scenarios_sunda2_java_segment_HS.csv',
          'random_scenarios_sunda2_sumatra_segment_HS.csv',
          'random_scenarios_sunda2_unsegmented_HS.csv')
        )
    # Names for each source representation in scenario_data [unsegmented, and
    # various segments]
    names_scenario_data = gsub('.csv', '', 
        gsub('random_scenarios_sunda2_', '', basename(scenario_data), 
             fixed=TRUE), 
        fixed=TRUE)

    # The multidomain directories for all SWALS model runs, for every scenario.
    # These do NOT need to be ordered in the same was as the scenario_data.
    tarred_multidomain_dirs = Sys.glob(paste0('../../swals/OUTPUTS/', run_series_name, 
        '/ptha18_random_scenarios_sunda2_row_*/RUN*'))

}else{
    stop(paste0('Unknown source_info: ', source_info))
}

# Output directory -- make sure it doesn't include "/"
output_directory = gsub('/', '-', run_series_name, fixed=TRUE)

# Compute the exceedance-rates for each of these depths (output in raster format)
# SEPARATELY FOR THE UNSEGMENTED SOURCE AND EACH OF THE SEGMENTS. The combined 
# source would be 0.5*(sum of unsegmented + union-of-segments), since in PTHA18
# they are weighted 50-50. 
depth_thresholds_for_exceedance_rate_calculations = 0.001  # c(0.001, 0.1, 0.5, 1, 3, 5, 8, 10)

# Given a source-model row index in the PTHA18 scenario database, find the
# SWALS tarred multidomain_dir that stores the tsunami model for that earthquake
# scenario. 
# This depends on the naming convention of the SWALS model output files, so 
# we make it a user input. Likely one will just need to edit the
# "matching_string" definition to conform to the model setup.
find_matching_md_dir<-function(row_indices, tarred_multidomain_dirs, source_zone){
    # Make a string with the start of the SWALS output folder name (beneath
    # ../../swals/OUTPUTS/...)
    matching_string = paste0('ptha18_random_scenarios_', source_zone, '_row_', 
        substring(as.character(1e+07 + row_indices), 2, 8), '_')

    # Match with the tarred_multidomain_dirs, with NA if we don't match or get multiple matches
    matching_ind = sapply(matching_string, f<-function(x){
        p = grep(x, tarred_multidomain_dirs)
        if(length(p) != 1) p = NA 
        return(p)})
    if(any(is.na(matching_ind))){
        stop('Could not find simulation matching scenario')
    }

    return(tarred_multidomain_dirs[matching_ind])
}

#
# END INPUTS
#

# Protect parallel calculations against failures of single tiles.
try_compute_exceedance_rates_on_tile<-function(
    input_domain_index_and_scenarios_name,
    tarred_multidomain_dirs, 
    scenario_databases, 
    output_directory, 
    depth_thresholds_for_exceedance_rate_calculations,
    raster_name_start){

    result = try(compute_exceedance_rates_on_tile(
        input_domain_index_and_scenarios_name,
        tarred_multidomain_dirs, 
        scenario_databases, 
        output_directory, 
        depth_thresholds_for_exceedance_rate_calculations,
        raster_name_start))
    return(list(result))
}

# Protect parallel calculations against failures of single tiles.
try_compute_exceedance_rates_and_error_variance_on_tile<-function(
    input_domain_index_and_scenarios_name,
    tarred_multidomain_dirs, 
    scenario_databases, 
    output_directory, 
    depth_thresholds_for_exceedance_rate_calculations,
    raster_name_start){

    result = try(compute_exceedance_rates_and_error_variance_on_tile(
        input_domain_index_and_scenarios_name,
        tarred_multidomain_dirs, 
        scenario_databases, 
        output_directory, 
        depth_thresholds_for_exceedance_rate_calculations,
        raster_name_start))
    return(list(result))
}




# Make space for the outputs
dir.create(output_directory, showWarnings=FALSE)

# Read the scenario databases
# Get all the csv data in a list with good names
scenario_databases = lapply(scenario_data, read.csv)
names(scenario_databases) = names_scenario_data

# For each scenario in the scenario_database, append the associated md_dir that
# holds the SWALS model outputs.
for(i in 1:length(scenario_databases)){
    scenario_databases[[i]]$md_dir = find_matching_md_dir(
        scenario_databases[[i]]$inds, tarred_multidomain_dirs, source_zone)
}

# Setup the parallel cluster. Note we will only do the exrate computation in
# parallel
local_cluster = makeCluster(MC_CORES)
clusterCall(local_cluster, fun=function(){library(utils); suppressPackageStartupMessages(library(rptha)) })
clusterExport(local_cluster, varlist=ls(all=TRUE))

# To run in parallel, make a list with each source representation and domain-index held together
sndi = expand.grid(names(scenario_databases), 1:max_domain_index)
parallel_job_list = vector(mode='list', length=nrow(sndi))
for(i in 1:nrow(sndi)){
    parallel_job_list[[i]] = list(
        domain_index = as.numeric(sndi[i,2]),
        scenarios_name = as.character(sndi[i,1]))
}
rm(sndi); gc()

# Run all the jobs
parallel_job_results = parLapplyLB(
    cl=local_cluster, 
    X=parallel_job_list, 
    #fun=try_compute_exceedance_rates_on_tile,
    fun=try_compute_exceedance_rates_and_error_variance_on_tile,
    tarred_multidomain_dirs=tarred_multidomain_dirs,
    scenario_databases=scenario_databases,
    output_directory=output_directory,
    depth_thresholds_for_exceedance_rate_calculations=depth_thresholds_for_exceedance_rate_calculations,
    raster_name_start=raster_name_start)

