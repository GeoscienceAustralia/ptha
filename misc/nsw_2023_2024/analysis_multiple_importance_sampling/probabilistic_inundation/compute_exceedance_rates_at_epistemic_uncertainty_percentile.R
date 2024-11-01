#
# Compute the exceedance-rate associated with a given depth threshold and epistemic-uncertainty percentile, as a raster.
#   Rscript compute_exceedance_rates_at_epistemic_uncertainty_percentile.R VARIABLE_NAME MD_BASE_DIRS DOMAIN_INDEX PERCENTILE_TO_USE EXCEEDANCE_THRESHOLD OUTPUT_DIR
# e.g.
#   Rscript compute_exceedance_rates_at_epistemic_uncertainty_percentile.R depth "../../swals/OUTPUTS/ptha18-NSW2023b-ID710.5-sealevel110cm/random_kermadectonga2/, ../../swals/OUTPUTS/ptha18-NSW2023-ID1315.5-sealevel110cm/random_kermadectonga2/, ../../swals/OUTPUTS/ptha18-NSW2023-ID4186.3-sealevel110cm/random_kermadectonga2/" 300 0.84 0.001 directory_to_hold_results
#
library(utils)
suppressPackageStartupMessages(library(rptha))

# Application specific information
asfm = new.env()
source('application_specific_file_metadata.R', local=asfm)

# Functions to help with epistemic uncertainty calculations
euf = new.env()
source('epistemic_uncertainty_functions.R', local=euf)

# Main ptha_access scripts
ptha18 = new.env()
source(asfm$get_ptha_results_script_file, local=ptha18, chdir=TRUE)

# R session resulting from "compute_rates_all_sources.R". This is
# needed to work with alternative logic-tree branches. It loads lots 
# of data so takes awhile to source.
ptha18_source_rate_env = new.env()
source(asfm$get_detailed_ptha18_source_zone_info_script_file,
       local=ptha18_source_rate_env, chdir=TRUE)

#browser()
#
# Key input parameters
#
inargs = commandArgs(trailingOnly=TRUE)
stopifnot(length(inargs) == 6)

# Folder containing a set of sub-folders of the form "ptha18*", which themselves contain raster_output_files.tar
VARIABLE_NAME = inargs[1] # 'depth' or 'max_stage'
MD_BASE_DIRS = gsub("^ ", "", scan(text=inargs[2], sep=",", what='character')) # Vector with multiple MD_BASE_DIRS for a single source-zone. Remove leading spaces if present.
DOMAIN_INDEX = as.numeric(inargs[3]) # 450
PERCENTILE_TO_USE = as.numeric(inargs[4]) # 0.84
EXCEEDANCE_THRESHOLD = as.numeric(inargs[5]) #0.001
OUTPUT_DIR = inargs[6] # 'my_output_directory_relative_to_here'

#
# Numerical parameters
#
# Optionally don't operate on every pixel (to speed calculations). E.g. a value
# of 3 will do proper calculations for the middle pixel in a 3x3 subgrid. The
# remaining pixels will later be filled using the middle pixel's value.
SUB_SAMPLE = 1
# How many random samples are used in the numerical percentile computation?
NRAND = 1e+04 
# Random seed used to make our 'random' function have deterministic results
# (required for root-finding). Try a different seed to check convergence
REPRODUCIBLE_SEED = 123  # 1234
# Put this value at sites (raster pixels) that should later be interpolated due to subsampling
# (distinct from 'missing data' which is a deliberate gap). Don't make it too 
# big, or a value that R will coerce to integer (and thus not be equal the 
# numeric equivalent).
NEEDS_INTERPOLATING = -9999.1

#
# End inputs
#

MC_CORES = asfm$DEFAULT_MC_CORES

# Inside the raster_output_files.tar archives, the rasters of interest should have names of the form
# paste0(raster_name_stub, DOMAIN_INDEX, '.tif')
raster_name_stub = asfm$get_raster_name_stub_from_variable_name(VARIABLE_NAME)

stopifnot(PERCENTILE_TO_USE > 0 & PERCENTILE_TO_USE < 1)
stopifnot(is.finite(c(EXCEEDANCE_THRESHOLD, DOMAIN_INDEX)))

dir.create(OUTPUT_DIR, showWarnings=FALSE)

# Get information on scenarios associated with this set of runs
sm = asfm$get_scenario_metadata_from_md_base_dir(MD_BASE_DIRS)
# Unpack the variables we use (comments show an example for sunda2)
source_info = sm$source_info # "random_sunda2"
source_zone = sm$source_zone # "sunda2"
scenario_base = sm$scenario_base # '../../sources/hazard/random_sunda2/'
raster_tar_files = sm$raster_tar_files # Sys.glob(paste0(scenario_base, 'ptha18*/raster_output_files.tar))

# Some metadata on the unsegmented/segments source representation
source_zone_segments_and_weights = 
    ptha18_source_rate_env$get_unsegmented_and_segmented_source_names_on_source_zone(source_zone)
# Names of the unsegmented source (first) and the segments, e.g.
#  c('sunda2', 'sunda2_arakan', 'sunda2_andaman', 'sunda2_sumatra', 'sunda2_java')
all_source_names = c(source_zone_segments_and_weights$unsegmented_name, 
    source_zone_segments_and_weights$segments) 
# Denote which indices in all_source_names are unsegmented vs segments
UNSEGMENTED_INDEX = 1
if(length(all_source_names) > 1){
    SEGMENTED_INDICES = seq(2, length(all_source_names))
}else{
    SEGMENTED_INDICES = c()
}
unsegmented_wt = source_zone_segments_and_weights$unsegmented_weight # 0.5
union_of_segments_wt = source_zone_segments_and_weights$union_of_segments_weight # 0.5

# List with only one entry named "logic_tree_mean_curve_HS" that contains a
# the pooled random scenarios, derived from a set of importance samples. In
# principle this doesn't need to be inside a list -- but it makes it easier to
# reuse some older code.
all_samples = sm$all_samples

# For the unsegmented/segments source representations, get some source frequency
# information for each logic tree branch, and other information used
# for subsequent calculations.
all_source_rate_info = list()
for(nm_i in all_source_names){
    all_source_rate_info[[nm_i]] = 
        # THIS IS NOT MODIFIED FOR MULTIPLE IMPORTANCE SAMPLING -- ADJUSTMENTS ARE FURTHER BELOW
        euf$get_key_source_information_all_logic_tree_branches_IS(
            source_zone, nm_i, 
            # The next line is the sampled scenarios in a data.frame --
            # identical for every source representation. We don't use the rates here directly.
            all_samples$logic_tree_mean_curve_HS, 
            ptha18, ptha18_source_rate_env)
}

# Convert the data to a structure suitable for pixel-by-pixel parallel calculation.
outputs = euf$make_all_pixel_data(all_samples[[1]], raster_tar_files, source_zone, 
    DOMAIN_INDEX, raster_name_stub, MC_CORES)
# Unpack the variables we use 
all_pixel_data = outputs$all_pixel_data # List with one entry per pixel, giving cell indices and vector with depths for raster_tar_files
scenarios_to_results_inds = outputs$scenarios_to_results_inds # Mapping between the scenario table and the raster results
output_matrix_dim = outputs$output_matrix_dim # Dimensions needed to store outputs
template_raster = outputs$template_raster # Template raster to store outputs
rm(outputs)
gc(verbose=FALSE)

# Compute the exrate for a hypothetical site (pixel) that is above the EXCEEDANCE_THRESHOLD for every
# scenario. This is a common case (e.g. in the ocean) and will be used in the parallel euf$get_exrate_percentile_at_pixel
# calculations to quickly get the solution (which improves the speed).
# Since we are doing multiple importance sampling, we require a different value for each sample
ALWAYS_WET_EXRATE = rep(NA, length=length(all_pixel_data[[1]]$sample_weights))
for(si in 1:length(ALWAYS_WET_EXRATE)){
    # Make some fake 'pixel' data 
    fake_weight = rep(0, length = length(ALWAYS_WET_EXRATE))
    fake_weight[si] = 1
    fake_wet_pixel_data = list(
        # Value is always > threshold
        model_runs_max_value = rep(EXCEEDANCE_THRESHOLD+1, length(all_samples[[1]]$inds)),
        # Pixel i/j indices ensure it is not skipped; counter unused here.
        i = floor(SUB_SAMPLE/2)+1, j = floor(SUB_SAMPLE/2)+1, counter=NA,
        sample_weights=fake_weight) 

    # Get the rate for an "always wet" pixel
    ALWAYS_WET_EXRATE[si] = euf$get_exrate_percentile_at_pixel(fake_wet_pixel_data,
        all_samples, all_source_rate_info, scenarios_to_results_inds,
        EXCEEDANCE_THRESHOLD, PERCENTILE_TO_USE,
        NRAND, SUB_SAMPLE, NEEDS_INTERPOLATING,
        NULL,  # Deliberate NULL where we'd usually pass ALWAYS_WET_EXRATE to trigger the calculation.
        UNSEGMENTED_INDEX, SEGMENTED_INDICES,
        unsegmented_wt, union_of_segments_wt,
        ptha18,
        REPRODUCIBLE_SEED)
}

# It is possible for ALWAYS_WET_EXRATE to be zero if we are considering
# a percentile that implies a zero event rate.
stopifnot(all(is.finite(ALWAYS_WET_EXRATE) & (ALWAYS_WET_EXRATE >= 0)))

# Fail-safe for parallel calculation of exceedance-rates for a single pixel
try_get_exrate_percentile_at_pixel<-function(pixel_data) try(euf$get_exrate_percentile_at_pixel(pixel_data,
    all_samples, all_source_rate_info, scenarios_to_results_inds,
    EXCEEDANCE_THRESHOLD, PERCENTILE_TO_USE,
    NRAND, SUB_SAMPLE, NEEDS_INTERPOLATING, ALWAYS_WET_EXRATE,
    UNSEGMENTED_INDEX, SEGMENTED_INDICES,
    unsegmented_wt, union_of_segments_wt,
    ptha18,
    REPRODUCIBLE_SEED))

#
# Note -- here if you just want to check a single point, you can use
# > single_point_value = try_get_exrate_percentile_at_pixel(all_pixel_data[[index_of_choice]])
# Likely index_of_choice corresponds to a point of interest. To find it you can get the cell
# indices i,j for every entry of all_pixel_data as:
# > i_s = unlist(lapply(all_pixel_data, function(x) x$i))
# > j_s = unlist(lapply(all_pixel_data, function(x) x$j))
# You can find these values at a point of interest with terra
# > tmp = extract(myrast, matrix(target_point, ncol=2), cells=TRUE)
# > j_ind = rowFromCells(myrast, tmp$cell)
# > i_ind = colFromCells(myrast, tmp$cell)
# > index_of_choice = which(i_s == i_ind & j_s == j_ind) ## Or I might have switched indices here
# > tmp[i_ind, j_ind] # Value of interest with 'terra' indexing (reverse order for 'stars')
#

#
# Clear unrequired memory before starting a parallel cluster
#
rm(ptha18_source_rate_env)
gc(verbose=FALSE)

# Setup a cluster
library(parallel)
local_cluster = makeCluster(MC_CORES)
ignore = clusterCall(local_cluster, fun=function(){ 
    library(utils); suppressPackageStartupMessages(library(rptha)) })
# Copy over all variables except 'all_pixel_data', which will be too large
vars_to_export = setdiff(ls(all=TRUE), 'all_pixel_data')
clusterExport(local_cluster, varlist=vars_to_export)

# Main calculation
all_pixel_results = parLapplyLB(cl = local_cluster, X = all_pixel_data, 
    fun = try_get_exrate_percentile_at_pixel, chunk.size=10)
stopCluster(local_cluster)

# Pack solution into a matrix
output_matrix = matrix(NA, ncol=output_matrix_dim[2], nrow=output_matrix_dim[1])
for(i in 1:length(all_pixel_results)){
    ind_i = all_pixel_data[[i]]$i
    ind_j = all_pixel_data[[i]]$j
    output_matrix[ind_i,ind_j] = all_pixel_results[[i]]
}

if(SUB_SAMPLE > 1){
    # Fill in the sub-sampled regions
    for(j in 1:ncol(output_matrix)){
        for(i in 1:nrow(output_matrix)){
            if(is.na(output_matrix[i,j])) next
            if(output_matrix[i,j] == NEEDS_INTERPOLATING){
                # Find the 'nearest' cell that was given values
                i_ind = floor((i-1)/SUB_SAMPLE)*SUB_SAMPLE + floor(SUB_SAMPLE/2) + 1
                j_ind = floor((j-1)/SUB_SAMPLE)*SUB_SAMPLE + floor(SUB_SAMPLE/2) + 1

                # Workaround for possible edge artefact, where cells would look to interpolate
                # from a cell that isn't there.
                if(i_ind > output_matrix_dim[1]) i_ind = i_ind - SUB_SAMPLE
                if(j_ind > output_matrix_dim[2]) j_ind = j_ind - SUB_SAMPLE

                output_matrix[i,j] = output_matrix[i_ind, j_ind]
            }
        }
    }
}

# Save to file
output_raster = setValues(template_raster, c(t(output_matrix)))
output_raster_filename = paste0(OUTPUT_DIR, '/', 
    source_info, '_', VARIABLE_NAME,
    '_rast_threshold_', EXCEEDANCE_THRESHOLD, 
    '_percentile_', 100*PERCENTILE_TO_USE, 
    '_subsam_', SUB_SAMPLE,
    '_Nrand_', NRAND,
    '_seed_', REPRODUCIBLE_SEED, 
    '_domain_index_', DOMAIN_INDEX, 
    '.tif')
writeRaster(output_raster, output_raster_filename, 
    options=c('COMPRESS=DEFLATE'), overwrite=TRUE)
