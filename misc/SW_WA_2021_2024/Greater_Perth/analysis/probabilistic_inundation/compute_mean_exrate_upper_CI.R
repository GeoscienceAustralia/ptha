
library(terra)

source_zones = c('outerrisesunda', 'sunda2')

##
## Initial case (later I revised the geometry)
##
#source_zone_dirs = paste0('ptha18-GreaterPerth-sealevel60cm-random_', source_zones, '-')
#summed_results_dir = 'ptha18-GreaterPerth-sealevel60cm-sum_of_source_zones'

#
# Run from inside reviseddomains_080422/highres_with_variance/
#
source_zone_dirs = paste0('ptha18-GreaterPerth-sealevel60cm-reviseddomain-highres-random_', source_zones, '-')
summed_results_dir = 'ptha18-GreaterPerth-sealevel60cm-reviseddomain-highres-sum_of_source_zones'
dir.create(summed_results_dir)

# Find the raster files and other relevant info
# @param sz_dir a source-zone dir, containing tifs to be used for calculations.
setup_files<-function(sz_dir){

    all_rasters = Sys.glob(paste0(sz_dir, '/*.tif'))

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
for(i in 1:length(source_zones)) sz_files[[source_zones[i]]] = setup_files(source_zone_dirs[i])

library(parallel)
mclapply(sz_files[[1]]$unique_domain_flags, sum_rasters_for_domain_i, sz_files=sz_files, output_dir=summed_results_dir,
    mc.cores=12)
