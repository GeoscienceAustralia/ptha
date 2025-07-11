# Compute the logic-tree-mean exceedance-rate and its Monte Carlo uncertainty
# by summing logic-tree-mean rasters over each source zone.
#
# NOTE: Previous versions of this script (for WA work) used separate rasters
# for the unsegmented and segmented source models on each source zone. In
# contrast, this version assumes that we have already computed the "combined"
# LTM solution on each source zone, so there is no need to combine
# unsegmented/segmented. The current version was developed for NSW work during
# my first usage of "importance sampling" (without any magnitude
# stratification) and fits with the outputs that are created in that context.
#
# Usage: Must run from directory containing "../application_specific_file_metadata.R"
#   Rscript compute_mean_exrate_upper_CI.R path/to/folder/with_subfolders_containing_results_on_each_source_zone/ variable threshold
# e.g. 
#   Rscript compute_mean_exrate_upper_CI.R ptha/sea_level_vary/highres_depth_with_variance/ depth 0.001
#   Rscript compute_mean_exrate_upper_CI.R ptha/sea_level_vary/highres_max_stage_with_variance/ max_stage 4.63
#
library(terra)

asfm = new.env()
source('../application_specific_file_metadata.R', local=asfm)

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

    # In this version of the code we should only be working with logic-tree-mean rasters
    stopifnot(all(startsWith(basename(all_rasters), "logic_tree_mean_HS_domain_")))

    unique_models = unique(raster_start)

    variance_rasters = all_rasters[grepl('__variance_of__exceedance_rate', all_rasters)]
    exrate_rasters = gsub('__variance_of__exceedance_rate', '__exceedance_rate', variance_rasters)
    stopifnot(all(file.exists(variance_rasters)))
    stopifnot(all(file.exists(exrate_rasters)))

    # How many domains?
    ND = length(exrate_rasters) #
    # Tag in filename that will identify the domain
    unique_domain_flags = paste0('_domain_', 1:ND, '_')

    # Double-check that each entry in unique_domain_flags match a single exrate raster. 
    count_matches = sapply(unique_domain_flags, function(x) sum(grepl(x, exrate_rasters)), USE.NAMES=FALSE)
    if(any(count_matches != 1)){
        print('Did not get a 1:1 match between unique_domain_flags and the LTM exceedance-rate rasters')
        print('This suggests a problem with identifying the rasters to operate on')
        stop('Deliberate halt')
    }
    rm(count_matches)

    return(environment())
}

#' Compute logic-tree-mean exceedance-rate by summing over the different source-zones. 
#' 
#' Our previous calculations have provided the logic-tree-mean solution for
#' individual source-zones. Here they are simply summed, taking care to treat
#' NA values well.
#'
#' @param unique_domain_flag_i a string like '_domain_i_' for the ith domain,
#' which will match the associated output tiff files. e.g. '_domain_24_' to
#' match the 24th domain.
#' @param sz_files The result of calling setup_files for the relevant
#' source-zone-directory
#' 
sum_rasters_for_domain_i<-function(unique_domain_flag_i, sz_files, output_dir){

    library(terra)
    # Loop over source-zones
    for(i in 1:length(sz_files)){

        sz = names(sz_files)[i]

        # Find rasters with the domain of interest in their name
        mtch = grep(unique_domain_flag_i, sz_files[[sz]]$exrate_rasters)
        mtch2 = grep(unique_domain_flag_i, sz_files[[sz]]$variance_rasters)

        # For the NSW case, length(mtch) should be 1. We have already computed
        # the logic-tree-mean result that implicitly includes the weighted sum
        # of unsegmented/segmented source representations.
        stopifnot(length(mtch) == 1)
        stopifnot(mtch2 == mtch)

        matching_exrate_rasters = sz_files[[sz]]$exrate_rasters[mtch]
        matching_variance_rasters = sz_files[[sz]]$variance_rasters[mtch]

        if(i == 1){
            # On the first loop iteration we make space to store the final results
            final_exrate = rast(matching_exrate_rasters[1])*0
            final_variance = rast(matching_variance_rasters[1])*0
            hasNA = is.na(final_exrate)
            final_exrate[hasNA] = 0
            final_variance[hasNA] = 0

            # Throw an error if there are any NA values
            stopifnot(!any(is.na(as.matrix(final_exrate))))
            stopifnot(!any(is.na(as.matrix(final_variance))))

            # Convenient to make output file names here
            exrate_outfile = paste0(output_dir, '/', 
                gsub('logic_tree_mean', 'summed', basename(matching_exrate_rasters[1])))
            variance_outfile = paste0(output_dir, '/', 
                gsub('logic_tree_mean', 'summed', basename(matching_variance_rasters[1])))
        }

        # Result on the source zone
        sz_exrate = rast(matching_exrate_rasters[1])
        sz_variance = rast(matching_variance_rasters[1])
        # Avoid problems with NA values by setting them to zero. The reason we
        # do that is that NA values denote sites that are "never wet" by
        # modelled events on the source zone (as well as halo regions from the
        # flow solver). But it is possible that some sites get wet from some
        # source zones, but not others. So if we sum results from multiple
        # source zones, we need to stop NA values from contaminating the
        # calculations. Later, we will reinstate NA values at sites that are not
        # affected by any source zone.
        hasNA = is.na(sz_exrate)
        sz_exrate[hasNA] = 0
        sz_variance[hasNA] = 0

        # Throw an error if there are any NA values
        stopifnot(!any(is.na(as.matrix(sz_exrate))))
        stopifnot(!any(is.na(as.matrix(sz_variance))))

        # Sum exceedance-rates over each source zone
        final_exrate = final_exrate + sz_exrate
        final_variance = final_variance + sz_variance

    }

    # Upper limit of 95% approximate confidence interval for the mean
    # exceedance-rate, accounting for the limited number of Monte Carlo
    # samples. Note this does NOT consider epistemic uncertainty
    exrate_upperCI = final_exrate + qnorm(0.975)*sqrt(final_variance)

    # Make zero regions NA, for ease of interpretation (and vrt creation) later on
    # This means NA regions are either "non_priority_domain regions" or else regions
    # that didn't get wet from any source zone.
    final_exrate[final_exrate == 0] = NA
    exrate_upperCI[exrate_upperCI == 0] = NA

    # It might be possible for the variance to be zero if the final_exrate is
    # non-zero (say if we only have 1 sample). So NA the final variance based
    # on the final exrate.
    final_variance[is.na(final_exrate)] = NA

    # Write files out
    writeRaster(final_exrate,   file=exrate_outfile, overwrite=TRUE, gdal=c('COMPRESS=DEFLATE'))
    writeRaster(final_variance, file=variance_outfile, overwrite=TRUE, gdal=c('COMPRESS=DEFLATE'))
    upper_CI_outfile = gsub('variance_of', 'Monte_Carlo_Upper_CI', variance_outfile)
    writeRaster(exrate_upperCI, file=upper_CI_outfile, overwrite=TRUE, gdal=c('COMPRESS=DEFLATE'))

    rm(final_exrate, final_variance, exrate_upperCI); gc()
        
    return(0)
}

create_vrt<-function(dir, threshold, variable=c("depth_as_max_stage_minus_elevation0")){
    system("module load gdal")
    ex_rate_cmd <- paste0("gdalbuildvrt -resolution highest ", 
        dir, "/exceedance_rate_with_threshold_", threshold, ".vrt ",
        dir, "/summed_HS_domain*",
        variable, 
        "_domain__exceedance_rate_with_threshold_", threshold, ".tif"
    )
    variance_cmd <- paste0("gdalbuildvrt -resolution highest ", 
        dir, "/variance_of__exceedance_rate_with_threshold_", threshold, ".vrt ",
        dir, "/summed_HS_domain*",
        variable,
        "_domain__variance_of__exceedance_rate_with_threshold_", threshold, ".tif"
    )
    upper_CI_cmd <- paste0("gdalbuildvrt -resolution highest ", 
        dir, "/Monte_Carlo_Upper_CI__exceedance_rate_with_threshold_", threshold, ".vrt ",
        dir, "/summed_HS_domain*",
        variable,
        "_domain__Monte_Carlo_Upper_CI__exceedance_rate_with_threshold_", threshold, ".tif"
    )
    system(ex_rate_cmd)
    system(variance_cmd)
    system(upper_CI_cmd)
}

# Get file info for the source zone
sz_files = list()
for(i in 1:length(source_zones)) sz_files[[source_zones[i]]] = setup_files(source_zone_dirs[i], threshold_of_interest)

library(parallel)
result = mclapply(sz_files[[1]]$unique_domain_flags, sum_rasters_for_domain_i, 
    sz_files=sz_files, output_dir=summed_results_dir,
    mc.cores=asfm$DEFAULT_MC_CORES)

create_vrt(summed_results_dir, threshold_of_interest, variable_of_interest)
