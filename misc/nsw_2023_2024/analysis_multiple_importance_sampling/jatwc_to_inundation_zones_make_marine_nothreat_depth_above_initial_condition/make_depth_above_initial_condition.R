#
# Compute the 'depth above the initial condition' raster corresponding to the
# largest "marine warning" and "no threat" scenarios. The zone is clipped to the
# corresponding polygon zones for consistency with what NSWSES have been using to date
#
# Use like (e.g. for Sydney-Coast marine_warning results)
#    Rscript make_depth_above_initial_conditions.R Sydney-Coast marine_warning
#

library(terra)

##
## INPUTS
##

input_parameters_group = commandArgs(trailingOnly=TRUE)
stopifnot(length(input_parameters_group) == 2)

# Which atws zone?
atws_zone = input_parameters_group[1] # e.g. "Sydney-Coast"
# Which warning type?
warning_type = input_parameters_group[2] # "marine_warning" or "no_threat"
stopifnot(any(warning_type %in% c('marine_warning', 'no_threat')))

# Polygon zone used to mask the outputs. Here we use the zones that were post
# processed, i.e. clipped to where the 1/2500 84% tsunami was below 1 cm,
# because this is what NSW SES have been using.
warning_zone_polygon_mask_shp = paste0('../jatwc_to_inundation_zones_edited/', atws_zone, '/',
    atws_zone, '_', warning_type, '_with-PTHA-exrate-limit_84pc_4e-04_where_1in2500_84pc_max_stage_exceeds_1.11/',
    atws_zone, '_', warning_type, '_with-PTHA-exrate-limit_84pc_4e-04_where_1in2500_84pc_max_stage_exceeds_1.11.shp')
stopifnot(file.exists(warning_zone_polygon_mask_shp))
warning_zone_polygon_mask = vect(warning_zone_polygon_mask_shp)

# Static sea level in model
MODEL_AMBIENT_SEALEVEL = 1.1 

# Skip sites with elevation below this.
IGNORE_SITES_WITH_ELEVATION_BELOW_M = 0.0

# In the final 'depth above initial condition', mask values below this. Tiny
# values can protect against floating point level roundoff (e.g. max-stage in a
# pond being 1e-08 above the initial condition due to numerical round off).
MASK_DEPTHS_BELOW_THIS_THRESHOLD = 0.0

# Find rasters with max stage in the warning zone
stage_rasts = Sys.glob(
    paste0('../jatwc_to_inundation/Inundation_zones/', atws_zone, '/', 
        warning_type, '_max_stage_domain*.tif'))
stopifnot(length(stage_rasts) > 0)
stopifnot(all(file.exists(stage_rasts)))

# Find the elevation rasters matching stage_rasts
matching_elev_rast_basenames = gsub(
    paste0(warning_type, "_max_stage"), "elevation0", basename(stage_rasts))
elevation_rasts = paste0('../../analysis_scenarios_ID710.5/jatwc_to_inundation/elevation_in_model/', 
    matching_elev_rast_basenames)
stopifnot(length(elevation_rasts) > 0)
stopifnot(all(file.exists(elevation_rasts)))

# Output directory
output_dir = paste0(atws_zone, '/', atws_zone, '_', warning_type, 
    '_depth_above_initial_condition_where_elevation_exceeds_0')

## END INPUTS

# output_dir must have a trailing '/'
if(!endsWith(output_dir, '/')) output_dir = paste0(output_dir, '/')

# Split into a list with one entry per raster, to make it easier to run in parallel
parallel_jobs = lapply(1:length(stage_rasts), 
    function(x) list(stage_rast_file=stage_rasts[x], elevation_rast_file=elevation_rasts[x]))

# Do the calculation of interest
compute_depth_with_care<-function(parallel_job){
    # stage raster
    sr = rast(parallel_job$stage_rast_file)
    # elevation raster
    er = rast(parallel_job$elevation_rast_file)

    # First estimate of the depth
    result = sr - er # Modified below

    # Remove sites outside the warning zone polygon mask
    result = mask(result, warning_zone_polygon_mask)

    # Remove sites with elevation below MSL (or whatever cutoff was specified)    
    sites_to_ignore = (er < IGNORE_SITES_WITH_ELEVATION_BELOW_M)
    result[sites_to_ignore] = NA

    # Region where we apply a depth adjustment (i.e. where depth above initial
    # condition is not the same as depth, and not ignored)
    sites_below_ambient_sl = (er < MODEL_AMBIENT_SEALEVEL) & (er >= IGNORE_SITES_WITH_ELEVATION_BELOW_M)

    # Depth above initial condition (will only be used at sites below ambient sea level)
    depth_above_initial_cond = (result - (MODEL_AMBIENT_SEALEVEL - er))
    depth_above_initial_cond = depth_above_initial_cond*(depth_above_initial_cond > 0) + 
        0*(depth_above_initial_cond <= 0) # Prevent negative depth

    result = 
        # This bit is the regular depth, masked to sites ABOVE ambient sea level
        result*(1-sites_below_ambient_sl) + 
        # This bit is the "depth above initial condition", masked to sites BELOW ambient sea level 
        depth_above_initial_cond*(sites_below_ambient_sl)

    # Remove newly introduced negative or zero depths, or depths below a given threshold
    negative = (result <= MASK_DEPTHS_BELOW_THIS_THRESHOLD)
    result[negative] = NA

    return(result)
}

# Do the calculations
result = lapply(parallel_jobs, compute_depth_with_care) # Serial
names(result) = paste0(atws_zone, '_', 
    gsub('max_stage', 
        'depth_above_initial_condition_where_elevation_exceeds_0',
        basename(stage_rasts), fixed=TRUE))

# Save the files
dir.create(output_dir, showWarnings=FALSE, recursive=TRUE)
for(i in 1:length(result)){
    output_rast = paste0(output_dir, names(result)[i])    
    writeRaster(result[[i]], output_rast, gdal=c('COMPRESS=DEFLATE'), overwrite=TRUE)
}

# Make a vrt
setwd(output_dir)
vrt_basename = paste0('all_', strsplit(names(result)[1], '_domain_')[[1]][1], '.vrt')
system(paste0('gdalbuildvrt -resolution highest ', vrt_basename, ' *.tif'))
