#
# Using the 'depth at 84 percentile' raster, estimate the 'max-stage at 84
# percentile' raster.
#
# This is cheaper than actually computing it, but NOT RECOMMENDED IN GENERAL.
# The estimate is limited because the depth rasters are CLIPPED to between
# 0-10m (for computational efficiency, since to calculate them we need root-finding).
#
# Note we might not have computed the depth for every domain tile
#
#
library(terra)


## INPUTS

# Rasters with the 84th percentile depth (clipped between 0 and 10m for calculation efficiency).
depth_at_84pc_rasts = Sys.glob('highres_domains_depth_at_epistemic_uncertainty_84pc/sum_of_all_sources_depth_rast_exrate_4e-04_percentile_84_subsam_1_Nrand_10000_seed_123_domain_index*.tif')
# Maximum depth that our epistemic uncertainty rasters report. Depths above
# this are "clipped" to CLIP_DEPTH.
CLIP_DEPTH = 9.99 # Almost 10 m (use slightly lower value, due to tolerance of root-finding used to make depth_at_84pc_rasts)
# In the outputs, put this value at sites where depth_at_84pc_rasts > CLIP_DEPTH
GREATER_THAN_CLIPDEPTH_FLAG = 1000

elevation_rasts = Sys.glob('../jatwc_to_inundation/elevation_in_model/elevation0_domain_*.tif')

output_dir = 'highres_domains_max_stage_clipped_if_depth_exceeds_10m_at_epistemic_uncertainty_84pc/'

## END INPUTS

# Get domain index for the elevation rasts (likely contains every model domain)
elevation_rasts_domain_index = sapply(elevation_rasts, function(x){
    s1 = strsplit(basename(x), '_domain_')[[1]][2]
    return(as.numeric(gsub(".tif", "", s1, fixed=TRUE)))
})
stopifnot(length(unique(elevation_rasts_domain_index)) == length(elevation_rasts_domain_index))   

# Get domain index for the depth@84 percentile rasts (likely doesn't contain
# every model domain, to reduce compute)
depth_rasts_domain_index = sapply(depth_at_84pc_rasts, function(x){
    s1 = strsplit(basename(x), '_domain_index_')[[1]][2]
    return(as.numeric(gsub(".tif", "", s1, fixed=TRUE)))
})
stopifnot(length(unique(depth_rasts_domain_index)) == length(depth_rasts_domain_index))

# Find matching indices
m1 = match(depth_rasts_domain_index, elevation_rasts_domain_index)
stopifnot(all(depth_rasts_domain_index == elevation_rasts_domain_index[m1]))

# Split into a list with one entry per raster, to make it easier to run in parallel
parallel_jobs = lapply(1:length(depth_at_84pc_rasts), 
    function(x) list(depth_rast_file=depth_at_84pc_rasts[x], elevation_rast_file=elevation_rasts[m1[x]]))

# Routine to compute depth, using a large numeric value to denote regions where the
# depth is > CLIP_DEPTH (because in these cases we just record a depth of 10m, so
# cannot reconstruct the stage).
sum_with_care<-function(parallel_job, greater_than_clipdepth_flag = GREATER_THAN_CLIPDEPTH_FLAG){
    dr = rast(parallel_job$depth_rast_file)
    er = rast(parallel_job$elevation_rast_file)

    result = dr + er    

    # The depth rasters never exceed 10m, so use a special flag to denote
    # regions with depth > 10 (which are invalid)
    deep_regions = (dr > CLIP_DEPTH)
    result = result * (1-deep_regions) + greater_than_clipdepth_flag*deep_regions

    # Make sure regions with zero depth show up as NA, rather than having
    # max-stage = elevation.
    zero_depth = (dr == 0)
    result[zero_depth] = NA

    return(result)
}

# Do the calculations
library(parallel)
#result = mclapply(parallel_jobs, sum_with_care, mc.cores=16) # Problems for some reason....?
result = lapply(parallel_jobs, sum_with_care) # Serial
names(result) = gsub('depth', 'max_stage_clipped_with_depth_exceeds_10m', basename(depth_at_84pc_rasts), fixed=TRUE)

# Save the files
dir.create(output_dir, showWarnings=FALSE)
for(i in 1:length(result)){
    output_rast = paste0(output_dir, names(result)[i])    
    writeRaster(result[[i]], output_rast, gdal=c('COMPRESS=DEFLATE'), overwrite=TRUE)
}
