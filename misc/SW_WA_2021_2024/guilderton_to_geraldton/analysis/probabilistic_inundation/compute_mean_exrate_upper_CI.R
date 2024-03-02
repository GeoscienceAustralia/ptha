# Compute the logic-tree-mean exceedance-rate and its Monte Carlo uncertainty
# by summing over existing rasters for the unsegmented/segmented models for
# each source zone.
#
# Usage: Must run from directory containing "application_specific_file_metadata.R"
#   Rscript compute_mean_exrate_upper_CI.R path/to/folder/with_subfolders_containing_results_on_each_source_zone/ variable threshold
# e.g. 
#   Rscript compute_mean_exrate_upper_CI.R ptha18-midwest-sealevel60cm/highres_with_variance depth 0.001
#
library(terra)

asfm = new.env()
source('application_specific_file_metadata.R', local=asfm)

# INPUT ARGUMENTS
folder_with_source_zone_logic_tree_mean_results = commandArgs(trailingOnly=TRUE)[1]
variable_of_interest = commandArgs(trailingOnly=TRUE)[2]
threshold_of_interest = commandArgs(trailingOnly=TRUE)[3]

source_zones = names(asfm$source_zone_modelled_tsunami_scenario_basedirs) # c('outerrisesunda', 'sunda2')

#
# Go to the folder containing logic-tree-mean results in sub-folders
# This might be a folder of the form
#   ptha18-BunburyBusseltonRevised-sealevel60cm/highres_with_variance
# created with make_directory_structure.sh
#
print(folder_with_source_zone_logic_tree_mean_results)
stopifnot(file.exists(folder_with_source_zone_logic_tree_mean_results))
setwd(folder_with_source_zone_logic_tree_mean_results)

matching_source_zone_dirs = lapply(source_zones, 
    function(x) Sys.glob(paste0('*-', variable_of_interest, '-LogicTreeMean-', x)))
source_zone_dirs = unlist(matching_source_zone_dirs)
stopifnot(all(unlist(lapply(matching_source_zone_dirs, length)) == 1)) # One directory per source-zone

# Make a directory to hold the summed results
summed_results_dir = gsub(source_zones[1], 'sum_of_source_zones', matching_source_zone_dirs[[1]])
dir.create(summed_results_dir, showWarnings=FALSE)

# Find the raster files and other relevant info
# @param sz_dir a source-zone dir, containing tifs to be used for calculations.
setup_files<-function(sz_dir, threshold_of_interest){

    all_rasters = Sys.glob(
        paste0(sz_dir, '/*exceedance_rate_with_threshold_', threshold_of_interest, '.tif'))

    raster_start = sapply(basename(all_rasters), 
        function(x) paste(strsplit(x, split="_")[[1]][1:2], collapse="_"), USE.NAMES=FALSE)

    unique_models = unique(raster_start)

    variance_rasters = all_rasters[grepl('__variance_of__exceedance_rate', all_rasters)]
    exrate_rasters = gsub('__variance_of__exceedance_rate', '__exceedance_rate', variance_rasters)
    stopifnot(all(file.exists(variance_rasters)))
    stopifnot(all(file.exists(exrate_rasters)))

    # How many domains? Count the unsegmented model (which occurs on all source-zones)
    ND = sum(grepl('unsegmented_HS_', exrate_rasters))
    # Tag in filename that will identify the domain
    unique_domain_flags = paste0('_domain_', 1:ND, '_')

    # Double-check that each entry in unique_domain_flags match a single exrate raster. 
    ui = grep('unsegmented_HS_', exrate_rasters)
    count_matches = sapply(unique_domain_flags, function(x) sum(grepl(x, exrate_rasters[ui])), USE.NAMES=FALSE)
    if(any(count_matches != 1)){
        print('Did not get a 1:1 match between unique_domain_flags and the unsegmented exceedance-rate rasters')
        print('This suggests a problem with identifying the rasters to operate on')
        stop('Deliberate halt')
    }
    rm(ui, count_matches)

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
        mtch = grep(unique_domain_flag_i, sz_files[[sz]]$exrate_rasters)
        matching_exrate_rasters = sz_files[[sz]]$exrate_rasters[mtch]
        matching_variance_rasters = sz_files[[sz]]$variance_rasters[mtch]

        if(length(mtch) > 1){
            # Site with mixed segmented/unsegmented treatment

            # In PTHA18
            # Result = 0.5 * unsegmented + 0.5 * (sum_of_segments)
            #        = 0.5 * sum_over_all_rasters
            local_exrate = rast(matching_exrate_rasters[1])*0
            local_variance = rast(matching_variance_rasters[1])*0

            # Avoid issues with NA values
            hasNA = is.na(local_exrate)
            local_exrate[hasNA] = 0
            local_variance[hasNA] = 0

            for(j in 1:length(mtch)){
                # Both variances and means are additive
                # But be careful with NA regions [the same region might have
                # inundation on another segment, need to account for that]

                # Mean exceedance-rate
                r1 = rast(matching_exrate_rasters[j])
                r1[is.na(r1)] = 0
                local_exrate = local_exrate + 0.5*r1

                # Variance of Monte-Carlo mean exceedance-rates
                r1 = rast(matching_variance_rasters[j])
                r1[is.na(r1)] = 0
                local_variance = local_variance + 0.5*r1
            }

        }else{
            # No segments -- unsegmented model has everything
            local_exrate = rast(matching_exrate_rasters[1])
            local_variance = rast(matching_variance_rasters[1])

            # Avoid issues with NA values
            hasNA = is.na(local_exrate)
            local_exrate[hasNA] = 0
            local_variance[hasNA] = 0

            # Convenient to make output file names here
            exrate_outfile = paste0(output_dir, '/', 
                gsub('unsegmented', 'summed', basename(matching_exrate_rasters[1])))
            # Convenient to make output file names here
            variance_outfile = paste0(output_dir, '/', 
                gsub('unsegmented', 'summed', basename(matching_variance_rasters[1])))

        }



        if(i == 1){
            # Make space to store summed results
            final_exrate = local_exrate
            final_variance = local_variance
        }else{
            # Sum results
            final_exrate = final_exrate + local_exrate
            final_variance = final_variance + local_variance
        }

    }

    rm(local_exrate, local_variance); gc()

    # Upper limit of 95% approximate confidence interval for the mean exceedance-rate, accounting for the
    # limited number of Monte Carlo samples. Note this does NOT consider epistemic uncertainty
    exrate_upperCI = final_exrate + qnorm(0.975)*sqrt(final_variance)

    # Make zero regions NA, for ease of interpretation later on
    final_exrate[final_exrate == 0] = NA
    #final_variance[final_variance == 0] = NA # Easier to keep variance zero
    exrate_upperCI[exrate_upperCI == 0] = NA

    writeRaster(final_exrate,   file=exrate_outfile, overwrite=TRUE, gdal=c('COMPRESS=DEFLATE'))
    writeRaster(final_variance, file=variance_outfile, overwrite=TRUE, gdal=c('COMPRESS=DEFLATE'))
    upper_CI_outfile = gsub('variance_of', 'Monte_Carlo_Upper_CI', variance_outfile)
    writeRaster(exrate_upperCI, file=upper_CI_outfile, overwrite=TRUE, gdal=c('COMPRESS=DEFLATE'))

    rm(final_exrate, final_variance, exrate_upperCI); gc()
        
    return(0)
}

# Get file info for the source zone
sz_files = list()
for(i in 1:length(source_zones)) sz_files[[source_zones[i]]] = setup_files(source_zone_dirs[i], threshold_of_interest)

library(parallel)
result = mclapply(sz_files[[1]]$unique_domain_flags, sum_rasters_for_domain_i, 
    sz_files=sz_files, output_dir=summed_results_dir,
    mc.cores=asfm$DEFAULT_MC_CORES)
