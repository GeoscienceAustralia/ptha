#
# After percentile exceedance-rate rasters have been computed, this script
# can sum them
# 
#     # 84th percentile calc
#     Rscript compute_sum_of_percentiles.R 84
#     # 16th percentile calc
#     Rscript compute_sum_of_percentiles.R 16
#

library(terra)

percentile_choice = commandArgs(trailingOnly=TRUE)[1]
stopifnot(any(percentile_choice == c('16', '84')))

if(percentile_choice == '84'){
    # 84th percentile
    source_zones = c('outerrisesunda', 'sunda2')
    source_zone_dirs = paste0('ptha18-BunburyBusseltonRevised-sealevel60cm/highres_epistemic_uncertainty/84pc/', 
        'ptha18-BunburyBusseltonRevised-sealevel60cm-depth_exrate_0.001_0.84_', source_zones, '/')
}else if(percentile_choice == '16'){
    # 16th percentile
    source_zones = c('outerrisesunda', 'sunda2')
    source_zone_dirs = paste0('ptha18-BunburyBusseltonRevised-sealevel60cm/highres_epistemic_uncertainty/16pc/', 
        'ptha18-BunburyBusseltonRevised-sealevel60cm-depth_exrate_0.001_0.16_', source_zones, '/')
}

summed_results_dir = gsub('outerrisesunda', 'sum_of_source_zones', source_zone_dirs[1])
dir.create(summed_results_dir, showWarnings=FALSE)

# Find the raster files and other relevant info
# @param sz_dir a source-zone dir, containing tifs to be used for calculations.
setup_files<-function(sz_dir){

    percentile_exrate_rasters = Sys.glob(paste0(sz_dir, '/*.tif'))

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
            output_raster_basename = gsub('outerrisesunda', 'sum_of_source_zones', basename(matching_exrate_rasters))
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
for(i in 1:length(source_zones)) sz_files[[source_zones[i]]] = setup_files(source_zone_dirs[i])

library(parallel)
mclapply(sz_files[[1]]$unique_domain_flags, sum_rasters_for_domain_i, sz_files=sz_files, output_dir=summed_results_dir,
    mc.cores=12)
