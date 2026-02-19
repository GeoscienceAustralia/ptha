#
# After percentile exceedance-rate rasters have been computed, this script
# can sum them
# 
#     # 84th percentile calc
#     Rscript compute_sum_of_percentiles.R ptha18-BunburyBusseltonRevised-sealevel60cm/highres_depth_epistemic_uncertainty 84 depth 0.001
#     # 16th percentile calc
#     Rscript compute_sum_of_percentiles.R ptha18-BunburyBusseltonRevised-sealevel60cm/highres_depth_epistemic_uncertainty 16 depth 0.001
#
library(terra)
asfm = new.env()
source('../application_specific_file_metadata.R', asfm)

# Directory containing one subdirectory for each percentile
results_folder = commandArgs(trailingOnly=TRUE)[1]

# The percentile we operate on (must match a percentile already created)
percentile_choice = commandArgs(trailingOnly=TRUE)[2]

# The variable and threshold we operate on
variable_of_interest = commandArgs(trailingOnly=TRUE)[3] # depth
threshold_of_interest = commandArgs(trailingOnly=TRUE)[4] # 0.001

source_zones = names(asfm$source_zone_modelled_tsunami_scenario_basedirs) # c('outerrisesunda', 'sunda2')

# For each source-zone, find the folder containing epistemic uncertainty outputs
szd_list = lapply(source_zones, function(sz){
    matching_folder = Sys.glob(paste0(results_folder, '/', percentile_choice, 'pc/*_', sz))
    if(length(matching_folder) != 1) stop(paste0('Did not find exactly 1 folder matching ', sz))
    return(matching_folder)
    })
source_zone_dirs = unlist(szd_list) # One directory per source zone

# Location for outputs
summed_results_dir = gsub(source_zones[1], 'sum_of_source_zones', source_zone_dirs[1])
dir.create(summed_results_dir, showWarnings=FALSE)

# Find the raster files and other relevant info
# @param sz_dir a source-zone dir, containing tifs to be used for calculations.
setup_files<-function(sz_dir, variable_of_interest, threshold_of_interest, percentile_choice){

    percentile_exrate_rasters = Sys.glob(
        paste0(sz_dir, '/*_', variable_of_interest, 
        '_rast_threshold_', threshold_of_interest, 
        '_percentile_', percentile_choice, 
        '_*.tif'))

    # How many domains? Count the unsegmented model (which occurs on all source-zones)
    ND = length(percentile_exrate_rasters)

    domain_IDS = sapply(percentile_exrate_rasters, FUN=function(x){
        fs = strsplit(x, '_')[[1]]
        lfs = length(fs)
        tag = gsub('.tif', '', fs[lfs], fixed=TRUE)
        return(tag)
        }, USE.NAMES=FALSE)

    # Tag in filename that will identify the domain
    unique_domain_flags = paste0('_domain_index_', domain_IDS, '.tif')

    return(environment())
}

#' Compute logic-tree-mean exceedance-rate by summing over the different source-zones. 
#' 
#' Within source-zones that have both unsegmented and segmented models, we
#' follow PTHA18 and have 50% weight on the unsegmented model and 50% weight
#' on the sub-of-segments (the calculation translates to just summing all the rasters and multiplying by 0.5).
#' For source-zones without segmented models we can use the unsegmented model alone to define the hazard.
#'
#' @param unique_domain_flag_i a string like '_domain_i_' for the ith domain, which will match the associated
#' output tiff files.
#' @param sz_files The result of calling setup_files for the relevant source-zone-directory
#' 
sum_rasters_for_domain_i<-function(unique_domain_flag_i, sz_files, output_dir){
    library(terra)

    # Loop over source-zones
    for(i in 1:length(sz_files)){

        sz = names(sz_files)[i]

        # Find rasters with the domain of interest in their name
        mtch = grep(unique_domain_flag_i, sz_files[[sz]]$percentile_exrate_rasters, fixed=TRUE)

        # Better not have more than one match
        stopifnot(length(mtch) <= 1)
        # We can have zero matches (e.g. did not do calculation for domain_1,
        # as it was large, and doesn't calculate inundation anyway).
        if(length(mtch) == 0){
            print(paste0('No match for ', unique_domain_flag_i))
            return(-1)
        }

        matching_exrate_rasters = sz_files[[sz]]$percentile_exrate_rasters[mtch]

        local_exrate = rast(matching_exrate_rasters[1])

        # Avoid issues with NA values
        hasNA = is.na(local_exrate)
        local_exrate[hasNA] = 0

        if(i == 1){
            # Make space to store summed results
            final_exrate = local_exrate

            # Convenient to make the output file basename here
            output_raster_basename = gsub(source_zones[1], 'sum_of_source_zones', basename(matching_exrate_rasters))
        }else{
            # Sum results
            final_exrate = final_exrate + local_exrate
        }

    }

    rm(local_exrate); gc()

    # Make zero regions NA, for ease of interpretation later on
    final_exrate[final_exrate == 0] = NA

    exrate_outfile = paste0(summed_results_dir, '/', output_raster_basename)
    writeRaster(final_exrate, file=exrate_outfile, overwrite=TRUE, gdal=c('COMPRESS=DEFLATE'))
    rm(final_exrate); gc()
        
    return(0)
}

# Get file info for the source zone
sz_files = list()
for(i in 1:length(source_zones)){
    sz_files[[source_zones[i]]] = setup_files(
        source_zone_dirs[i], variable_of_interest, threshold_of_interest, percentile_choice)
}

# Sanity check
if(length(source_zones) > 1){
    for(i in 2:length(source_zones)){
        if(length(sz_files[[i]]$percentile_exrate_rasters) != length(sz_files[[1]]$percentile_exrate_rasters)){
            stop('Error in creating or finding rasters: There is a different number of percentile rasters on each source')
        }
    }
}
library(parallel)
result = mclapply(sz_files[[1]]$unique_domain_flags, sum_rasters_for_domain_i, 
    sz_files=sz_files, output_dir=summed_results_dir,
    mc.cores=asfm$DEFAULT_MC_CORES)
