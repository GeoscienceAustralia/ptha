#
# The max-stage @ 1/2500 84th percentile raster has an undesirable property,
# due to details of our "variable agnostic" calculations of the threshold at a given
# exceedance-rate and epistemic uncertainty.
#
# The issue is
# - At sites that are dry at the chosen return period, the code returns the
#   lower-bound of the solution search space, which is typically equal to the
#   background sea level (1.1m in our case).
# - This value could be below the bed elevation!
#
# It is desirable to cut these areas from the output rasters to prevent any misinterpretation.
# 
# Notice this isn't as much of a a problem for "depth" because the lower bound of the search space is 0.0,
# which clearly means "dry". But it is still useful to set depth=0.0 to depth=NA, to make it
# easier to correctly visualise the results.
library(terra)

## INPUTS
max_stage_at_84pc_rasts = Sys.glob(
    'highres_domains_max_stage_at_epistemic_uncertainty_84pc/sum_of_all_sources_max_stage_rast_exrate_4e-04_percentile_84_subsam_1_Nrand_10000_seed_123_domain_index*.tif')

# Lower limit of seach space (plus slight tolerance to avoid floating-point equality issues). Typically the background sea level.
remove_below_this_value = 1.10001

elevation_rasts = Sys.glob('../jatwc_to_inundation/elevation_in_model/elevation0_domain_*.tif')

output_dir = 'highres_domains_max_stage_at_epistemic_uncertainty_84pc_masked_in_dry_areas/'
## END INPUTS


# Get domain index for the elevation rasts (likely contains every model domain)
elevation_rasts_domain_index = sapply(elevation_rasts, function(x){
    s1 = strsplit(basename(x), '_domain_')[[1]][2]
    return(as.numeric(gsub(".tif", "", s1, fixed=TRUE)))
})
stopifnot(length(unique(elevation_rasts_domain_index)) == length(elevation_rasts_domain_index))   

# Get domain index for the max_stage@84 percentile rasts (likely doesn't contain
# every model domain, to reduce compute)
max_stage_rasts_domain_index = sapply(max_stage_at_84pc_rasts, function(x){
    s1 = strsplit(basename(x), '_domain_index_')[[1]][2]
    return(as.numeric(gsub(".tif", "", s1, fixed=TRUE)))
})
stopifnot(length(unique(max_stage_rasts_domain_index)) == length(max_stage_rasts_domain_index))

# Find matching indices
m1 = match(max_stage_rasts_domain_index, elevation_rasts_domain_index)
stopifnot(all(max_stage_rasts_domain_index == elevation_rasts_domain_index[m1]))

# Split into a list with one entry per raster, to make it easier to run in parallel
parallel_jobs = lapply(1:length(max_stage_at_84pc_rasts), 
    function(x) list(max_stage_rast_file=max_stage_at_84pc_rasts[x], elevation_rast_file=elevation_rasts[m1[x]]))

clip_max_stage_dry_regions<-function(parallel_job){
    # depth raster
    sr = rast(parallel_job$max_stage_rast_file)
    # elevation raster
    er = rast(parallel_job$elevation_rast_file)

    # Remove sites matching the lower-limit
    lower_limit_region = (sr < (0*sr + remove_below_this_value))
    sr[lower_limit_region] = NA

    #dry_depth = (sr <= er)
    #sr[dry_depth] = NA

    return(sr)
}

result = lapply(parallel_jobs, clip_max_stage_dry_regions)
names(result) = gsub('max_stage', 
    'max_stage_masked_in_dry_regions',
    basename(max_stage_at_84pc_rasts), fixed=TRUE)

# Save the files
dir.create(output_dir, showWarnings=FALSE)
for(i in 1:length(result)){
    output_rast = paste0(output_dir, names(result)[i])    
    writeRaster(result[[i]], output_rast, gdal=c('COMPRESS=DEFLATE'), overwrite=TRUE)
}
