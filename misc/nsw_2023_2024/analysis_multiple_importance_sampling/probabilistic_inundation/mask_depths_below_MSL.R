# Compute "masked" depth rasters that only show the depth at sites above some elevation.
# Run with
#    Rscript mask_depths_below_MSL.R 1in2500_84pc
# e.g.
#    Rscript "nsw_full_coast_MIS_highres_domains_depth_at_epistemic_uncertainty_84pc/*_rast_exrate_4e-04_*.tif" "../../analysis_scenarios_ID710.5/jatwc_to_inundation/elevation_in_model/*.tif" nsw_full_coast_MIS_highres_domains_depth_where_elevation_exceeds_0_84pc_4e-04
#
# Be sure to use "quotes" to prevent the wildcards being expanded on the shell.
# Also be sure to match only the tif files of interest (e.g. in case there are
# rasters at multiple exceedance-rates in the one folder).
#
# This lets us restrict the depth to onshore sites.
#
# The script also sets sites with depth=0 to depth=NA, to prevent an "easy to
# make mistake" in data visualisation, although nowadays that is taken care of
# in an earlier post-processing step using 
#     tidy_lower_bounds_in_thresholds_at_exceedance_rate_of_epistemic_uncertainty_percentile.R
#

library(terra)

## INPUTS
input_par = commandArgs(trailingOnly=TRUE)[1]
if(input_par == "") input_par = '1in2500_84pc'

if(input_par == '1in2500_84pc'){
    string_matching_depth_rasters = 'nsw_full_coast_MIS_max_depth_1in2500_84pc/*.tif'
    string_matching_elevation_rasters = "../../analysis_scenarios_ID710.5/jatwc_to_inundation/elevation_in_model/*.tif"
    output_dir = 'nsw_full_coast_MIS_max_depth_masked_at_elevation_below_0_1in2500_84pc/'
}else if(input_par == '1in250_50pc'){
    string_matching_depth_rasters = 'nsw_full_coast_MIS_max_depth_1in250_50pc/*.tif'
    string_matching_elevation_rasters = "../../analysis_scenarios_ID710.5/jatwc_to_inundation/elevation_in_model/*.tif"
    output_dir = 'nsw_full_coast_MIS_max_depth_masked_at_elevation_below_0_1in250_50pc/'
}else{
    stop(paste0('unknown input_par ', input_par))
}

IGNORE_SITES_WITH_ELEVATION_BELOW_M = 0.0 # Skip sites with elevation below this.
MASK_DEPTHS_BELOW_THRESHOLD = -10 # No longer needed since I use the 'tidy....' step previously #0.003 # Mask depths less than this depth {reflecting accuracy limits in the uniroot estimate of the depth}

# Rasters with the 84th percentile depth (clipped between 0 and 10m for
# calculation efficiency).
depth_at_84pc_rasts = Sys.glob(string_matching_depth_rasters)
#    'highres_domains_depth_at_epistemic_uncertainty_84pc/sum_of_all_sources_depth_rast_exrate_4e-04_percentile_84_subsam_1_Nrand_10000_seed_123_domain_index*.tif')

elevation_rasts = Sys.glob(string_matching_elevation_rasters) #'../jatwc_to_inundation/elevation_in_model/elevation0_domain_*.tif')

#output_dir = 'highres_domains_depth_at_epistemic_uncertainty_84pc_masked_at_elevation_below_zero/'
if(!endsWith(output_dir, '/')) output_dir = paste0(output_dir, '/')

## END INPUTS

stopifnot(length(depth_at_84pc_rasts) > 0)
stopifnot(length(elevation_rasts) > 0)

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
    dr[dr < MASK_DEPTHS_BELOW_THRESHOLD] = NA

    # Get rid of results where the elevation is too low
    to_remove = (er < IGNORE_SITES_WITH_ELEVATION_BELOW_M)
    dr[to_remove] = NA

    return(dr)

}

result = lapply(parallel_jobs, compute_masked_depth)
names(result) = gsub('depth', 
    'depth_where_elevation_exceeds_zero',
    basename(depth_at_84pc_rasts), fixed=TRUE)

# Save the files
dir.create(output_dir, showWarnings=FALSE)
for(i in 1:length(result)){
    output_rast = paste0(output_dir, names(result)[i])    
    writeRaster(result[[i]], output_rast, gdal=c('COMPRESS=DEFLATE'), overwrite=TRUE)
}
