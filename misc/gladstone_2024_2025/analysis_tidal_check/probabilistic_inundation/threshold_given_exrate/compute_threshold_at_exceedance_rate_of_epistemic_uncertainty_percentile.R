#
# Compute the depth associated with a given exceedance-rate and epistemic-uncertainty percentile, as a raster.
# Run like:
#   Rscript compute_threshold_at_exceedance_rate_of_epistemic_uncertainty_percentile.R VARIABLE_NAME DOMAIN_INDEX PERCENTILE_TO_USE EXCEEDANCE_RATE MINIMUM_DEPTH MAXIMUM_DEPTH DEPTH_TOL OUTPUT_DIR
# e.g.:
#   Rscript compute_threshold_at_exceedance_rate_of_epistemic_uncertainty_percentile.R depth 450 0.84 0.0004 0.005 10.0 0.004 directory_to_hold_results
# or (max_stage with a data-defined lower bound for the uniroot search)
#   Rscript compute_threshold_at_exceedance_rate_of_epistemic_uncertainty_percentile.R max_stage 450 0.84 0.0004 adaptive_minimum 10.6 0.001 directory_to_hold_results
# or (max_speed with a data-defined upper bound for the uniroot search)
#   Rscript compute_threshold_at_exceedance_rate_of_epistemic_uncertainty_percentile.R max_speed 450 0.84 0.0004 0.0 adaptive_maximum 0.001 directory_to_hold_results
#
library(utils)
suppressPackageStartupMessages(library(rptha))

# Application specific information
asfm = new.env()
source('../application_specific_file_metadata.R', local=asfm)

# Functions to help with epistemic uncertainty calculations
euf = new.env()
source('../epistemic_uncertainty_functions.R', local=euf)

# Main ptha_access scripts
ptha18 = new.env()
source(asfm$get_ptha_results_script_file, local=ptha18, chdir=TRUE)

# R session resulting from "compute_rates_all_sources.R". This is
# needed to work with alternative logic-tree branches. It loads lots 
# of data so takes awhile to source.
ptha18_source_rate_env = new.env()
source(asfm$get_detailed_ptha18_source_zone_info_script_file,
       local=ptha18_source_rate_env, chdir=TRUE)


#
# Key input parameters
#
inargs = commandArgs(trailingOnly=TRUE)
stopifnot(length(inargs) == 8)

# Folder containing a set of sub-folders of the form "ptha18*", which themselves contain raster_output_files.tar
VARIABLE_NAME = inargs[1] # 'depth' or 'max_stage'
#MD_BASE_DIR = inargs[2] #'../../../swals/OUTPUTS/ptha18-GreaterPerth2023-sealevel60cm/random_outerrisesunda/'
DOMAIN_INDEX = as.numeric(inargs[2]) # 450
PERCENTILE_TO_USE = as.numeric(inargs[3]) # 0.84
EXCEEDANCE_RATE = as.numeric(inargs[4]) # 1/2500
if(inargs[5] == 'adaptive_minimum'){
    # Infer MINIMUM DEPTH search range from the data. In this case, if the
    # smallest non-NA tsunami in a given pixel corresponds to an exceedance-rate that is rarer than
    # EXCEEDANCE_RATE, then returned pixel value will be NA. This is typically what we want for
    # areas that are dry at the specified exceedance-rate.
    MINIMUM_DEPTH = 'adaptive_minimum'
}else{
    MINIMUM_DEPTH = as.numeric(inargs[5]) # 0.0005
}
if(inargs[6] == 'adaptive_maximum'){
    # Infer MAXIMUM DEPTH search range from the data.
    MAXIMUM_DEPTH = 'adaptive_maximum'
}else{
    MAXIMUM_DEPTH = as.numeric(inargs[6]) # 10.0
}
DEPTH_TOL = as.numeric(inargs[7]) # 0.005 - Tolerance for the depth used in the root-finding function
OUTPUT_DIR = inargs[8] # 'my_output_directory_relative_to_here'

#
# Numerical parameters
#
# Optionally don't operate on every pixel (to speed calculations). E.g. a value
# of 3 will do proper calculations for the middle pixel in a 3x3 subgrid. The
# remaining pixels will later be filled using the middle pixel's value.
SUB_SAMPLE = 1 #3
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
stopifnot(all(is.finite(c(EXCEEDANCE_RATE, DOMAIN_INDEX, DEPTH_TOL))))
stopifnot(all(!is.na(c(MINIMUM_DEPTH, MAXIMUM_DEPTH))))

dir.create(OUTPUT_DIR, showWarnings=FALSE)

prepare_source_zone_data<-function(MD_BASE_DIR){
    # On each source zone, set up the data required for exceedance-rate calculations with epistemic uncertainty,
    # and make functions that do that.

    # Get information on scenarios associated with this set of runs
    sm = asfm$get_scenario_metadata_from_md_base_dir(MD_BASE_DIR)
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
    # filename with the random scenarios, sampled with importance sampling. In
    # principle this doesn't need to be inside a list -- but it makes it easier to
    # reuse some older code. FIXME: Better naming here and in the config code.
    all_source_samples = sm$all_source_samples 

    # Read the random scenarios. Here length(all_source_names)==1 so this needn't
    # be inside a list, but putting it in the list helps me reuse older code. 
    stopifnot(length(all_source_samples) == 1)
    all_samples = list(source_zone = read.csv(all_source_samples[[1]]))
    names(all_samples) = names(all_source_samples)[1]
    # Ensure the number of raster tar files matches the number of unique sampled scenarios
    stopifnot(length(raster_tar_files) == length(unique(all_samples[[1]]$inds)))

    # For the unsegmented/segments source representations, get some source frequency
    # information for each logic tree branch, and other information used
    # for subsequent calculations.
    all_source_rate_info = list()
    for(nm_i in all_source_names){
        all_source_rate_info[[nm_i]] = 
            euf$get_key_source_information_all_logic_tree_branches_IS(
                source_zone, nm_i, 
                # The next line is the sampled scenarios in a data.frame --
                # identical for every source representation. We don't use the rates here directly.
                all_samples$logic_tree_mean_curve_HS, 
                ptha18, ptha18_source_rate_env)
    }

    # Convert the data to a structure suitable for pixel-by-pixel parallel calculation.
    outputs = euf$make_all_pixel_data(
        all_samples[[1]]$inds,
        raster_tar_files,
        source_zone, 
        DOMAIN_INDEX,
        raster_name_stub,
        MC_CORES,
        MINIMUM_DEPTH,
        MAXIMUM_DEPTH
    )
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
    #
    # Make some fake 'pixel' data
    fake_threshold = 10
    fake_wet_pixel_data = list(
        # Value is always > threshold
        model_runs_max_value = rep(fake_threshold+1, length(all_samples[[1]]$inds)),
        # Pixel i/j indices ensure it is not skipped; counter unused here.
        i = floor(SUB_SAMPLE/2)+1, j = floor(SUB_SAMPLE/2)+1, counter=NA) 
    # Get the rate for an "always wet" pixel
    ALWAYS_WET_EXRATE = euf$get_exrate_percentile_at_pixel(fake_wet_pixel_data,
        all_samples, all_source_rate_info, scenarios_to_results_inds,
        fake_threshold, PERCENTILE_TO_USE,
        NRAND, SUB_SAMPLE, NEEDS_INTERPOLATING,
        NULL,  # Deliberate NULL where we'd usually pass ALWAYS_WET_EXRATE to trigger the calculation.
        UNSEGMENTED_INDEX, SEGMENTED_INDICES,
        unsegmented_wt, union_of_segments_wt,
        ptha18,
        REPRODUCIBLE_SEED)
    # It is possible for ALWAYS_WET_EXRATE to be zero if we are considering
    # a percentile that implies a zero event rate.
    stopifnot(is.finite(ALWAYS_WET_EXRATE) & (ALWAYS_WET_EXRATE >= 0))

    # Fail-safe for parallel calculation of exceedance-rates for a single pixel
    try_get_exrate_percentile_at_pixel<-function(pixel_data, chosen_threshold) try(
        euf$get_exrate_percentile_at_pixel(pixel_data,
            all_samples, all_source_rate_info, scenarios_to_results_inds,
            chosen_threshold, PERCENTILE_TO_USE,
            NRAND, SUB_SAMPLE, NEEDS_INTERPOLATING, ALWAYS_WET_EXRATE,
            UNSEGMENTED_INDEX, SEGMENTED_INDICES,
            unsegmented_wt, union_of_segments_wt,
            ptha18,
            REPRODUCIBLE_SEED))

    return(environment())
}

# For each sourcezone, prepare the data/functions
source_zone_data = lapply(asfm$source_zone_modelled_tsunami_scenario_basedirs,
    function(x) prepare_source_zone_data(x))

# Remove "all_pixel_data" from each entry of source_zone_data (and move it to the current core)
# This prevents memory blow-outs in parallel. We will distribute the data over all cores.
all_source_zone_pixel_data = lapply(source_zone_data, function(x) x$all_pixel_data)
for(i in 1:length(source_zone_data)){
    source_zone_data[[i]]$all_pixel_data = NULL
}
# Check that all_source_zone_pixel_data has the same number of cells on each source zone. 
# It should, because the entries correspond to cells on the raster corresponding to DOMAIN_INDEX.
all_lens = unlist(lapply(all_source_zone_pixel_data, length))
stopifnot(all(all_lens == all_lens[1]))

# Restructure "all_source_zone_pixel_data" to make parallel distribution cleaner
# We want a list with 1 entry per pixel.
# Each entry contains data for ALL source zones
tmp_source_zone_pixel_data = vector(mode='list', length=all_lens[1])
for(i in 1:length(tmp_source_zone_pixel_data)){
    tmp_source_zone_pixel_data[[i]] = lapply(all_source_zone_pixel_data, function(x) x[[i]])    
}
all_source_zone_pixel_data = tmp_source_zone_pixel_data
rm(tmp_source_zone_pixel_data); gc()

#' A function to mimise for root-finding.
#'
#' This assumes comonotonic dependence between the uncertainties on each
#' source zone (i.e. a conservative model of dependence). 
#'
#' @param depth a value of the threshold
#' @param source_zone_pixel_data a single entry of all_source_zone_pixel_data
#' (which has one entry for each source-zone)
#' @return The difference between the summed exceedance-rate and EXCEEDANCE_RATE.
#' 
rootfun<-function(depth, source_zone_pixel_data){

    stopifnot(length(source_zone_pixel_data) == length(source_zone_data) &
              all(names(source_zone_pixel_data) == names(source_zone_data)))

    # Get the exceedance-rate with the chosen depth/percentile on each source-zone separately
    exrates = rep(NA, length=length(source_zone_data))
    for(sz_i in 1:length(source_zone_data)){
        exrates[sz_i] = source_zone_data[[sz_i]]$try_get_exrate_percentile_at_pixel(
            source_zone_pixel_data[[sz_i]], 
            depth)
    }

    if(!any(is.na(exrates))){
        if(all(exrates == NEEDS_INTERPOLATING)){
            # Case that the pixel should be interpolated from neighbours
            return(NEEDS_INTERPOLATING)
        }
    }else{
        # There is at least one NA exrate
        #
        # It is possible that one source zone has exrates[i] = NA (i.e. dry all the time)   
        # but another does not. In that case, we should treat the exceedance-rate on the 'dry'
        # source zone as zero.
        if(any(is.finite(exrates))) exrates[is.na(exrates)] = 0
    }


    # Comontonic assumption to justify the sum(exrates) here.
    result = sum(exrates) - EXCEEDANCE_RATE 

    return(result)
}

# Function to find the exceedance-rate of interest via root finding.
find_threshold_at_exrate<-function(source_zone_pixel_data){

    # Minimum search depth over all source zones. Nontrivial if MINIMUM_DEPTH=='adaptive_minimum'
    min_d = min(unlist(lapply(source_zone_pixel_data, function(x) x$min))) - DEPTH_TOL/3
    # Maximum search depth over all source zones. Nontrivial if MAXIMUM_DEPTH=='adaptive_maximum'
    max_d = max(unlist(lapply(source_zone_pixel_data, function(x) x$max))) + DEPTH_TOL/3

    f_lower = rootfun(min_d, source_zone_pixel_data)
    f_upper = rootfun(max_d, source_zone_pixel_data)

    if(is.na(f_lower) & is.na(f_upper)) return(NA)

    # At least one of f_lower, f_upper is not NA.
    # Then neither should be NA
    if(is.na(f_lower) | is.na(f_upper)) stop('This should not happen!')

    # Special case where a pixel should be interpolated from neighbours later on.
    if(f_lower == NEEDS_INTERPOLATING){
        stopifnot(f_upper == NEEDS_INTERPOLATING)
        return(NEEDS_INTERPOLATING)
    }

    # We expect (f_lower > 0) and (f_upper < 0). Otherwise the
    # value of interest is outside our bounds, so quick exit
    if(f_lower < 0){
        if(MINIMUM_DEPTH == 'adaptive_minimum'){
            # Most likely the cell is "dry" at this exceedance rate. Or we 
            # requested an exceedance-rate such that our smallest scenario
            # has an (individual) occurrence-rate that is larger than the exceedance-rate.
            return(NA) 
        }else{
            return(min_d)
        }
    }
    if(f_upper > 0) return(max_d) # The real maxima is higher, clipped

    # Find the root
    findroot = uniroot(rootfun, 
        source_zone_pixel_data = source_zone_pixel_data, 
        lower=min_d, upper=max_d, tol=DEPTH_TOL)

    return(findroot$root)
}

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
vars_to_export = setdiff(ls(all=TRUE), 'all_source_zone_pixel_data')
clusterExport(local_cluster, varlist=vars_to_export)

# Main calculation
all_pixel_results = parLapplyLB(
    cl = local_cluster,
    X = all_source_zone_pixel_data, 
    fun = find_threshold_at_exrate,
    chunk.size=50
)
stopCluster(local_cluster)

# Pack solution into a matrix
output_matrix_dim = source_zone_data[[1]]$output_matrix_dim
output_matrix = matrix(NA, ncol=output_matrix_dim[2], nrow=output_matrix_dim[1])
for(i in 1:length(all_pixel_results)){
    ind_i = all_source_zone_pixel_data[[i]][[1]]$i
    ind_j = all_source_zone_pixel_data[[i]][[1]]$j
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
template_raster = source_zone_data[[1]]$template_raster
output_raster = setValues(template_raster, c(t(output_matrix)))
output_raster_filename = paste0(OUTPUT_DIR, '/', 
    'sum_of_all_sources_', VARIABLE_NAME,
    '_rast_exrate_', EXCEEDANCE_RATE, 
    '_percentile_', 100*PERCENTILE_TO_USE, 
    '_subsam_', SUB_SAMPLE,
    '_Nrand_', NRAND,
    '_seed_', REPRODUCIBLE_SEED, 
    '_domain_index_', DOMAIN_INDEX, 
    '.tif')
writeRaster(output_raster, output_raster_filename, 
    options=c('COMPRESS=DEFLATE'), overwrite=TRUE)
