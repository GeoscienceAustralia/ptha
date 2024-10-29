#
#
#
library(terra)
library(parallel)

#
# INPUTS
#

source_zones = c("alaskaaleutians",
    "kermadectonga2",
    "newhebrides2",
    "outerrise_kermadectonga",
    "outerrise_puysegur",
    "outerrisenewhebrides",
    "puysegur2",
    "solomon2",
    "southamerica")

# We use a single weight raster for each source zone, and its name must match this tag
# Surround by underscores to prevent matches of e.g. '1000' to '10000'
weight_raster_exrate_tag =  "_1in10000_" 

# Output folder -- we will append the source_zone to this.
# Make sure it ends with `-VARIABLE-LogicTreeMean-`, for instance '-depth-LogicTree-Mean-'
output_folder_startname = './ptha18-NSW2023-MIS-sealevel110cm/highres_depth_with_variance/ptha18-NSW2023-MIS-sealevel110cm-depth-LogicTreeMean-'

MC_CORES = 16 # How many cores
skip_domain_1 = TRUE # Should we skip domain 1 (global)? Often much faster and less memory hungry to avoid the global domain.

# Get the mean rate rasters for a source zone
get_mean_raster_files<-function(source_zone){
    # Rates for each tile
    mean_raster_files = list(
        'scenarios_ID710.5' = 
            Sys.glob(paste0('../../analysis_scenarios_ID710.5/probabilistic_inundation/ptha18-NSW2023b-ID710.5-sealevel110cm/highres_depth_with_variance/ptha18-NSW2023b-ID710.5-sealevel110cm-depth-LogicTreeMean-', 
                source_zone, '/*_domain__exceedance_rate_with_threshold_0.001*.tif')),
        'scenarios_ID1315.5' = 
            Sys.glob(paste0('../../analysis_scenarios_ID1315.5/probabilistic_inundation/ptha18-NSW2023-ID1315.5-sealevel110cm/highres_depth_with_variance/ptha18-NSW2023-ID1315.5-sealevel110cm-depth-LogicTreeMean-', 
                source_zone, '/*_domain__exceedance_rate_with_threshold_0.001*.tif')),
        'scenarios_ID4186.3' = 
            Sys.glob(paste0('../../analysis_scenarios_ID4186.3/probabilistic_inundation/ptha18-NSW2023-ID4186.3-sealevel110cm/highres_depth_with_variance/ptha18-NSW2023-ID4186.3-sealevel110cm-depth-LogicTreeMean-', 
                source_zone, '/*_domain__exceedance_rate_with_threshold_0.001*.tif'))
        )
    return(mean_raster_files)
}

# Get the variance rasters for a source zone
get_variance_raster_files<-function(source_zone){
    # Variances for each tile

    mean_raster_files = get_mean_raster_files(source_zone)
    variance_raster_files = lapply(mean_raster_files, function(x){ 
        gsub("_domain__exceedance_rate_with_threshold_", 
            "_domain__variance_of__exceedance_rate_with_threshold_", x, fixed=TRUE)}
        )
    return(variance_raster_files)
}

# Get the importance sample weight rasters for each source zone
get_weight_rasters<-function(source_zone, weight_raster_exrate_tag){
    # Spatially variable weights for each importance sample.
    # One raster per importance sample.

    mean_raster_files = get_mean_raster_files(source_zone)
    weight_raster_files = lapply(names(mean_raster_files), function(x){
        rast(paste0('../../sources/hazard/multiple_importance_sampling_weights/importance_sample_weight_rasters/weight_raster_', 
            source_zone, '_r', weight_raster_exrate_tag, x, '.tif'))})
    names(weight_raster_files) = names(mean_raster_files)
    return(weight_raster_files)
}

#
# END INPUTS
#

weight_rasters_sum_to_1<-function(weight_raster_files){
    # Weight rasters must sum to 1
    w1 = weight_raster_files[[1]]
    for(i in 2:length(weight_raster_files)) w1 = w1 + weight_raster_files[[i]]
    w1_vals = as.matrix(w1)
    result = (all(abs(w1_vals - 1) < 1e-06))
    return(result) 
}

extract_domain_index<-function(raster_filename){
    # Get domain index from input raster filename
    as.numeric(strsplit(strsplit(basename(raster_filename), "logic_tree_mean_HS_domain_")[[1]][2],  
        "_depth_as_max_stage_minus")[[1]][1])
} 

# Workhorse function to do calculations for a single domain
weight_rasters<-function(raster_files_mean, raster_files_var, weight_raster_files, output_file_mean, output_file_var, output_file_95CI, skip_domain_1){

    sample_names = names(raster_files_mean)
    stopifnot(names(raster_files_mean) == sample_names)
    stopifnot(names(raster_files_var) == sample_names)
    stopifnot(names(weight_raster_files) == sample_names)

    # If we should skip domain 1, then do it
    if(skip_domain_1 & grepl("_domain_1_", raster_files_mean[1])) return(invisible(0))

    # Get weight rasters with the same geometry as the other files
    wts = lapply(as.list(raster_files_mean), function(x) rast(x)*0)
    default_weight = 1/length(sample_names)
    for(sn in sample_names){
        wts[[sn]] = resample( weight_raster_files[[sn]], wts[[sn]], method='bilinear')
        # In NA regions, default to equal weights for all samples
        wts[[sn]][is.na(wts[[sn]])] = default_weight
    }

    # Sum the rasters. Start with a template blank raster
    r_mean = rast(raster_files_mean[[1]]) * 0
    r_mean[is.na(r_mean)] = 0 # Remove NA's because they just represent 'not wet', but different samples might have non-NA values there
    r_var = r_mean
    r_mean_upper95CI = r_mean
    for(i in 1:length(raster_files_mean)){

        # Accumulate weighted mean, sum_i { mean_i * w_i }
        r_i = rast(raster_files_mean[[i]])
        r_i[is.na(r_i)] = 0 # Remove NA's
        r_i = r_i * wts[[i]]
        r_mean = r_mean + r_i 

        # Accumulate variance, sum_i { mean_i * w_i**2 }
        r_i = rast(raster_files_var[[i]])
        r_i[is.na(r_i)] = 0 # Remove NA's as above
        r_i = r_i * wts[[i]] * wts[[i]] # Squared weights
        r_var = r_var + r_i

    }
    # 95% confidence interval upper limit, assuming a normal distribution
    r_mean_upper95CI = r_mean + qnorm(0.975) * sqrt(r_var)

    # Reintroduce NA values
    r_mean[r_mean == 0] = NA
    r_var[r_var == 0] = NA
    r_mean_upper95CI[r_mean_upper95CI == 0] = NA

    writeRaster(r_mean, filename = output_file_mean, gdal=c('COMPRESS=DEFLATE'), overwrite=TRUE)
    writeRaster(r_var, filename = output_file_var, gdal=c('COMPRESS=DEFLATE'), overwrite=TRUE)
    writeRaster(r_mean_upper95CI, filename = output_file_95CI, gdal=c('COMPRESS=DEFLATE'), overwrite=TRUE)

    rm(list=ls()); gc()
    
    return(invisible(0))
}

for(source_zone in source_zones){
    print(source_zone)

    output_folder = paste0(output_folder_startname, source_zone, '/')
    dir.create(output_folder, recursive=TRUE, showWarnings=FALSE)

    mean_raster_files = get_mean_raster_files(source_zone)
    variance_raster_files = get_variance_raster_files(source_zone)
    weight_raster_files = get_weight_rasters(source_zone, weight_raster_exrate_tag)

    # Names for output rasters
    output_raster_files_mean = paste0(output_folder, '/', basename(mean_raster_files[[1]]))
    output_raster_files_variance = paste0(output_folder, '/', basename(variance_raster_files[[1]]))
    # For the upper 95% confidence interval, the naming convention follows that used for individual samples.
    output_raster_files_95CI = gsub(
        '_domain__variance_of__exceedance_rate_with_threshold_', 
        '_domain__Monte_Carlo_Upper_CI__exceedance_rate_with_threshold_', 
        output_raster_files_variance, fixed=TRUE)

    # Check that weight rasters sum to 1
    if(!weight_rasters_sum_to_1(weight_raster_files)){
        stop('Weight rasters do not sum to 1')
    }

    # Check domain indices are aligned in all the input files
    check_domain_indices_align<-function(){
        # Get domain indices for all input files, and check they align
        domain_indices_mean = lapply(mean_raster_files, function(x) 
            unlist(lapply(x, extract_domain_index)))
        domain_indices_var = lapply(variance_raster_files, function(x) 
            unlist(lapply(x, extract_domain_index)))
        for(i in 2:length(domain_indices_mean)){
            stopifnot(all(domain_indices_mean[[i]] == domain_indices_mean[[1]]))
            stopifnot(all(domain_indices_var[[i]] == domain_indices_mean[[1]]))
        }
    }
    check_domain_indices_align()

    # Convert key inputs to data.frames, so we can extract along rows
    mean_raster_files = as.data.frame(mean_raster_files)
    variance_raster_files = as.data.frame(variance_raster_files)

    parallel_fun<-function(i){
        try(weight_rasters(mean_raster_files[i,], variance_raster_files[i,], weight_raster_files,
            output_raster_files_mean[i], output_raster_files_variance[i], output_raster_files_95CI[i], skip_domain_1))
    }

    parallel_job = mclapply(1:nrow(mean_raster_files), parallel_fun, mc.cores=MC_CORES, mc.preschedule=FALSE)
    is_try_error = unlist(lapply(parallel_job, function(x) is(x, 'try-error')))
    if(any(is_try_error)){
        print(parallel_job)
    }
}
