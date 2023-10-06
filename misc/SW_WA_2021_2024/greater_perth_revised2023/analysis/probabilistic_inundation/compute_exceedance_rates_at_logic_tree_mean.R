#
# Compute rasters with rate (depth > threshold).
# Can do similar computations for max-stage, arrival-time, etc.
#
library(utils)
suppressPackageStartupMessages(library(rptha))
suppressPackageStartupMessages(library(parallel))
asfm = new.env()
source('application_specific_file_metadata.R', local=asfm)
errc = new.env()
source('exceedance_rate_raster_calculations.R', local=errc)

# Print this if input arguments are problematic
usage_msg = paste0(
    'Usage example \n',
    '  Using runs in "../../swals/OUTPUTS/ptha18_tonga_MSL0/", for domains 1:4, to compute rates of (depth > 0.001) and separately (depth > 1.2), the command is:\n',
    '  Rscript compute_exceedance_rates_for_threshold_depth_logic_tree_mean_newParallelPartition.R ../../swals/OUTPUTS/ptha18_tonga_MSL0/ 4 depth 0.001 1.2')

#
# INPUTS
#

input_args = commandArgs(trailingOnly=TRUE)
if(!(length(input_args) >= 4)) stop(usage_msg)

# The name of the folder containing many tsunami model runs
MD_BASE_DIR = input_args[1] 

# The maximum domain index (corresponding to an integer in the tif filename).
# We do calculations for domain 1, 2, 3, ... MAX_DOMAIN_INDEX
MAX_DOMAIN_INDEX = as.numeric(input_args[2]) 

# Either 'depth' or 'max_stage' or 'arrival_time' or ...
VARIABLE_NAME = input_args[3] 

# Compute the exceedance-rates for each of these thresholds (output in raster
# format) SEPARATELY FOR THE UNSEGMENTED SOURCE AND EACH OF THE SEGMENTS.
EXCEEDANCE_THRESHOLDS = as.numeric(input_args[-c(1,2,3)])

# Sanity check
if(!(all(EXCEEDANCE_THRESHOLDS > 0) & (MAX_DOMAIN_INDEX > 0))) stop(usage_msg)

# Number of cores -- the computation is distributed over these.
MC_CORES = asfm$DEFAULT_MC_CORES
if(MC_CORES > detectCores()){
    stop(' MC_CORES exceeds the number of cores on your machine. Bad idea.')
}

# Rasters that we operate on must have filename starting with this.
raster_name_start = asfm$get_raster_name_stub_from_variable_name(VARIABLE_NAME)

# Get information on scenarios associated with this set of runs
sm = asfm$get_scenario_metadata_from_md_base_dir(MD_BASE_DIR)
# Unpack the variables we use (comments show an example for sunda2)
source_info = sm$source_info # "random_sunda2"
source_zone = sm$source_zone # "sunda2"
scenario_base = sm$scenario_base # '../../sources/hazard/random_sunda2/'
#all_source_names = sm$all_source_names # c('sunda2', 'sunda2_arakan', 'sunda2_andaman', 'sunda2_sumatra', 'sunda2_java')
all_source_samples = sm$all_source_samples # List with on entry per all_source_names, containing the random scenarios
#UNSEGMENTED_INDEX = sm$UNSEGMENTED_INDEX # 1 
#SEGMENTED_INDICES = sm$SEGMENTED_INDICES # c(2,3,4,5)
#unsegmented_wt = sm$unsegmented_wt # 0.5
#union_of_segments_wt = sm$segmented_wt # 0.5
#raster_tar_files = sm$raster_tar_files # Sys.glob(paste0(MD_BASE_DIR, '/ptha18*/raster_output_files.tar))
tarred_multidomain_dirs = sm$tarred_multidomain_dirs # Sys.glob(paste0(MD_BASE_DIR, '/ptha18*/RUN*.tar))

# Files with scenario row indices and other metadata
scenario_data = unlist(all_source_samples, use.names=FALSE)

# Names for each source representation in scenario_data [unsegmented, and
# various segments] This is different to all_source_names. There is no strong
# reason but easier not to change.
names_scenario_data = gsub('.csv', '',
    gsub(paste0('random_scenarios_', source_zone, '_'), '', basename(scenario_data),
        fixed=TRUE),
    fixed=TRUE)

# Output directory -- make sure it doesn't include "/"
output_directory = paste0(basename(dirname(MD_BASE_DIR)), '-', VARIABLE_NAME, '-LogicTreeMean-', source_zone) # "ptha18-GreaterPerth2023-sealevel60cm-depth-LogicTreeMean-outerrisesunda"

#
# END INPUTS
#

# Protect parallel calculations against failures of single tiles.
try_compute_exceedance_rates_and_error_variance_on_tile<-function(
    input_domain_index_and_scenarios_name_and_exrate,
    tarred_multidomain_dirs, 
    scenario_databases, 
    output_directory, 
    #EXCEEDANCE_THRESHOLDS, # Now contained in input_domain_index_and_scenarios_name_and_exrate
    raster_name_start){

    result = try(errc$compute_exceedance_rates_and_error_variance_on_tile(
        input_domain_index_and_scenarios_name_and_exrate,
        tarred_multidomain_dirs, 
        scenario_databases, 
        output_directory, 
        input_domain_index_and_scenarios_name_and_exrate$exceedance_threshold, #EXCEEDANCE_THRESHOLDS,
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
    scenario_databases[[i]]$md_dir = asfm$find_matching_md_data(
        scenario_databases[[i]]$inds, tarred_multidomain_dirs, source_zone)
}

# Setup the parallel cluster. Note we will only do the exrate computation in
# parallel
local_cluster = makeCluster(MC_CORES)
clusterCall(local_cluster, fun=function(){library(utils); suppressPackageStartupMessages(library(rptha)) })
clusterExport(local_cluster, varlist=ls(all=TRUE))

# To run in parallel, make a list with each source representation and domain-index and exceedance-threshold held together
sndi = expand.grid(names(scenario_databases), 1:MAX_DOMAIN_INDEX, EXCEEDANCE_THRESHOLDS)
parallel_job_list = vector(mode='list', length=nrow(sndi))
for(i in 1:nrow(sndi)){
    parallel_job_list[[i]] = list(
        domain_index = as.numeric(sndi[i,2]),
        scenarios_name = as.character(sndi[i,1]),
        exceedance_threshold = as.numeric(sndi[i,3]))
}
rm(sndi); gc()

# Run all the jobs
parallel_job_results = parLapplyLB(
    cl=local_cluster, 
    X=parallel_job_list, 
    fun=try_compute_exceedance_rates_and_error_variance_on_tile,
    tarred_multidomain_dirs=tarred_multidomain_dirs,
    scenario_databases=scenario_databases,
    output_directory=output_directory,
    # EXCEEDANCE_THRESHOLDS=EXCEEDANCE_THRESHOLDS, # Now contained in parallel_job_list
    raster_name_start=raster_name_start)
stopCluster(local_cluster)
