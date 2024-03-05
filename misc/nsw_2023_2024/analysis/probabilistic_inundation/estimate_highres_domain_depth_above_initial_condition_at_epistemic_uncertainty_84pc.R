# Using the 'depth at 84 percentile' raster, estimate the 'depth above the initial condition'
# raster.
#
# The initial condition for the NSW model is 1.1m. Some areas are flooded even
# at this depth (e.g. inside estuaries where the tides are attenuated). At such sites
# the tsunami is often small (due to attenuation) and the raw depth values can 
# greatly overestimate the size of the tsunami. For example
#     Elevation = 0.5m 
#     Initial_condition = 1.1m
#     Max_stage = 1.15m
#     So the tsunami is 5cm high, but the raw depth will be reported as 0.65.
# For this reason it is of interest to consider the "depth above the initial condition" 
# (5cm in this example). Basically this reduces the "artificially high" depths in
# sites flooded by the initial condition, without affecting depths at other sites.
#
# Note we might not have computed the depth for every domain tile
#
library(terra)

## INPUTS

MODEL_AMBIENT_SEALEVEL = 1.1 # Static sea level in model

IGNORE_SITES_WITH_ELEVATION_BELOW_M = 0.0 # Skip sites with elevation below this.

# Rasters with the 84th percentile depth (clipped between 0 and 10m for calculation efficiency).
depth_at_84pc_rasts = Sys.glob('highres_domains_depth_at_epistemic_uncertainty_84pc/sum_of_all_sources_depth_rast_exrate_4e-04_percentile_84_subsam_1_Nrand_10000_seed_123_domain_index*.tif')
# Maximum depth that the depth_at_84pc_rasts report. Greater depths were 
# "clipped" when those rasters were created
CLIP_DEPTH = 9.99 # Almost 10 m (slightly less due to tolerance of root-finding used to make depth_at_84pc_rasts)

elevation_rasts = Sys.glob('../jatwc_to_inundation/elevation_in_model/elevation0_domain_*.tif')

output_dir = 'highres_domains_depth_above_initial_condition_at_epistemic_uncertainty_84pc/'

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

# Do the calculation of interest
sum_with_care<-function(parallel_job, greater_than_10m_flag = 1000){
    # Depth raster
    dr = rast(parallel_job$depth_rast_file)
    # elevation raster
    er = rast(parallel_job$elevation_rast_file)

    result = dr # Modified eblow

    # Make sure regions with zero depth show up as NA, rather than having
    # max-stage = elevation.
    zero_depth = (dr == 0)
    result[zero_depth] = NA

    # Remove sites with elevation below MSL (or whatever cutoff was specified)    
    sites_to_ignore = (er < IGNORE_SITES_WITH_ELEVATION_BELOW_M)
    result[sites_to_ignore] = NA

    # Region where we apply a depth adjustment
    sites_below_ambient_sl = (er < MODEL_AMBIENT_SEALEVEL) & (er >= IGNORE_SITES_WITH_ELEVATION_BELOW_M)

    adjusted_depth = (result - (MODEL_AMBIENT_SEALEVEL - er))
    adjusted_depth = adjusted_depth*(adjusted_depth > 0) + 0*(adjusted_depth <= 0)

    result = result*(1-sites_below_ambient_sl) + 
        # This bit subtracts the 'non-tsunami-depth' from sites below the ambient sea level.
        adjusted_depth*(sites_below_ambient_sl)

    return(result)
}

# Do the calculations
library(parallel)
#result = mclapply(parallel_jobs, sum_with_care, mc.cores=16) # Problems for some reason....?
result = lapply(parallel_jobs, sum_with_care) # Serial
names(result) = gsub('depth', 'depth_above_initial_condition_where_elevation_exceeds_0', basename(depth_at_84pc_rasts), fixed=TRUE)

# Save the files
dir.create(output_dir, showWarnings=FALSE)
for(i in 1:length(result)){
    output_rast = paste0(output_dir, names(result)[i])    
    writeRaster(result[[i]], output_rast, gdal=c('COMPRESS=DEFLATE'), overwrite=TRUE)
}
