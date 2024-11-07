# Compute the 'depth above the initial condition' raster, using the
# - depth @ 84 percentile rasters
# - elevation rasters

#
# Why?
# The initial condition for the NSW model is 1.1m. Some areas are flooded by
# the initial condition (e.g. inside estuaries where the tides are attenuated).
# At such sites the tsunami is often small (due to attenuation) and the raw
# depth values can greatly overestimate the size of the tsunami. For example
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

input_parameters_group = commandArgs(trailingOnly=TRUE)
if(input_parameters_group == "") input_parameters_group = '1in2500_84pc'

if(input_parameters_group == '1in2500_84pc'){
    ## INPUTS
    string_matching_depth_rasters = "nsw_full_coast_MIS_max_depth_1in2500_84pc/*_rast_exrate_4e-04_*.tif" 
    string_matching_elevation_rasters = "../../analysis_scenarios_ID710.5/jatwc_to_inundation/elevation_in_model/*.tif" 
    output_dir = 'nsw_full_coast_MIS_max_depth_above_initial_condition_1in2500_84pc_at_sites_with_elevation_above_0'
}else if(input_parameters_group == "1in250_50pc"){
    # Inputs for another case
    string_matching_depth_rasters = "nsw_full_coast_MIS_max_depth_1in250_50pc/*_rast_exrate_0.004_*.tif" 
    string_matching_elevation_rasters = "../../analysis_scenarios_ID710.5/jatwc_to_inundation/elevation_in_model/*.tif" 
    output_dir = 'nsw_full_coast_MIS_max_depth_above_initial_condition_1in250_50pc_at_sites_with_elevation_above_0'
}else{
    stop(paste0('unknown input parameters group: ', input_parameters_group))
}

#string_matching_depth_rasters = commandArgs(trailingOnly=TRUE)[1]
#string_matching_elevation_rasters = commandArgs(trailingOnly=TRUE)[2]
#output_dir = commandArgs(trailingOnly=TRUE)[3]

MODEL_AMBIENT_SEALEVEL = 1.1 # Static sea level in model

IGNORE_SITES_WITH_ELEVATION_BELOW_M = 0.0 # Skip sites with elevation below this.

MASK_DEPTHS_BELOW_THRESHOLD = -10 # No longer needed since the tidy_....R script takes care of this # Mask depths less than this depth {reflecting accuracy limits in the uniroot estimate of the depth}

# Rasters with 84th percentile max-depth.
depth_at_84pc_rasts = Sys.glob(string_matching_depth_rasters)

elevation_rasts = Sys.glob(string_matching_elevation_rasters)

#output_dir = # 'highres_domains_depth_above_initial_condition_at_epistemic_uncertainty_84pc/'

## END INPUTS

stopifnot(length(depth_at_84pc_rasts) > 0)
stopifnot(length(elevation_rasts) > 0)
if(!endsWith(output_dir, '/')) output_dir = paste0(output_dir, '/')

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
compute_depth_with_care<-function(parallel_job){
    # depth raster
    dr = rast(parallel_job$depth_rast_file)
    # elevation raster
    er = rast(parallel_job$elevation_rast_file)

    # First estimate of the depth
    result = dr # Modified below

    # Remove negative or zero depths (or depths within this range from the uniroot tolerance)
    negative = (result < MASK_DEPTHS_BELOW_THRESHOLD)
    result[negative] = NA

    # Remove sites with elevation below MSL (or whatever cutoff was specified)    
    sites_to_ignore = (er < IGNORE_SITES_WITH_ELEVATION_BELOW_M)
    result[sites_to_ignore] = NA

    # Region where we apply a depth adjustment
    sites_below_ambient_sl = (er < MODEL_AMBIENT_SEALEVEL) & (er >= IGNORE_SITES_WITH_ELEVATION_BELOW_M)

    # Depth above initial condition (will only be used at sites below ambient
    # sea level)
    depth_above_initial_cond = (result - (MODEL_AMBIENT_SEALEVEL - er))
    depth_above_initial_cond = depth_above_initial_cond*(depth_above_initial_cond > 0) + 
        0*(depth_above_initial_cond <= 0) # Prevent negative depth

    result = 
        # This bit is the regular depth, masked to sites ABOVE ambient sea level
        result*(1-sites_below_ambient_sl) + 
        # This bit is the "depth above initial condition", masked to sites BELOW ambient sea level 
        depth_above_initial_cond*(sites_below_ambient_sl)

    # Remove newly introduced negative or zero depths. Deliberately avoid MASK_DEPTHS_BELOW_THRESHOLD here,
    # because we already accounted for it above.
    negative = (result <= 0.0)
    result[negative] = NA

    return(result)
}

# Do the calculations
library(parallel)
#result = mclapply(parallel_jobs, compute_depth_with_care, mc.cores=16) # Problems for some reason....?
result = lapply(parallel_jobs, compute_depth_with_care) # Serial
names(result) = gsub('depth', 
    'depth_above_initial_condition_where_elevation_exceeds_0',
    basename(depth_at_84pc_rasts), fixed=TRUE)

# Save the files
dir.create(output_dir, showWarnings=FALSE)
for(i in 1:length(result)){
    output_rast = paste0(output_dir, names(result)[i])    
    writeRaster(result[[i]], output_rast, gdal=c('COMPRESS=DEFLATE'), overwrite=TRUE)
}
