#
# Compute "masked" depth rasters that only show the depth at sites above some elevation.
# 
# This lets us restrict the depth to onshore sites.
#

library(terra)

## INPUTS

IGNORE_SITES_WITH_ELEVATION_BELOW_M = 0.0 # Skip sites with elevation below this.

# Rasters with the 84th percentile depth (clipped between 0 and 10m for
# calculation efficiency).
depth_at_84pc_rasts = Sys.glob(
    'highres_domains_depth_at_epistemic_uncertainty_84pc/sum_of_all_sources_depth_rast_exrate_4e-04_percentile_84_subsam_1_Nrand_10000_seed_123_domain_index*.tif')

elevation_rasts = Sys.glob('../jatwc_to_inundation/elevation_in_model/elevation0_domain_*.tif')

output_dir = 'highres_domains_depth_at_epistemic_uncertainty_84pc_masked_at_elevation_below_zero/'

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

compute_masked_depth<-function(parallel_job){
    # depth raster
    dr = rast(parallel_job$depth_rast_file)
    # elevation raster
    er = rast(parallel_job$elevation_rast_file)

    # Get rid of zero depths
    dr[dr == 0] = NA

    # Get rid of results where the elevation is too low
    to_remove = (er < IGNORE_SITES_WITH_ELEVATION_BELOW_M)
    dr[to_remove] = NA

    return(dr)

}

result = lapply(parallel_jobs, compute_masked_depth)
names(result) = gsub('depth', 
    'depth_above_initial_condition_where_elevation_exceeds_0',
    basename(depth_at_84pc_rasts), fixed=TRUE)

# Save the files
dir.create(output_dir, showWarnings=FALSE)
for(i in 1:length(result)){
    output_rast = paste0(output_dir, names(result)[i])    
    writeRaster(result[[i]], output_rast, gdal=c('COMPRESS=DEFLATE'), overwrite=TRUE)
}
