#
# Compute the maximum modelled depth over all of the 'marine-warning' scenarios
# combined.
#
# This can only be run AFTER the "map_threat_levels_in_zone.R" script has been run
#

#library(stars)
library(terra)

# Get a set of rasters covering the zone -- here we choose no_threat, although any choice would be fine.
# The purpose is just to easily identify which domains should be included.
ATWS_ZONE_NAME = commandArgs(trailingOnly=TRUE)[1]
ATWS_Zone_name_nospace = gsub(' ', '-', ATWS_ZONE_NAME)

# Find tif files with the max-stage (max over all marine warning scenarios). Note
# these tifs are NOT limited by the PTHA, although that is unlikely to have much effect.
all_marine_warning_tifs = Sys.glob(paste0('Inundation_zones/', ATWS_Zone_name_nospace, 
    '/marine_warning_max_stage*.tif'))
stopifnot(all(file.exists(all_marine_warning_tifs)))

# Elevation files for all the tifs
all_elevation_tifs = paste0('./elevation_in_model/', 
    gsub('marine_warning_max_stage_', 'elevation0_', basename(all_marine_warning_tifs)))
stopifnot(all(file.exists(all_elevation_tifs)))

# Write tifs here
output_dir = paste0(dirname(all_marine_warning_tifs[1]), '/', ATWS_Zone_name_nospace, 
    '_marine_warning_depth_maxima')
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)

# Compute the max-depth
process_file_i<-function(i){

    #max_stage_union = read_stars(all_marine_warning_tifs[i])
    max_stage_union = rast(all_marine_warning_tifs[i])
    #elevation = read_stars(all_elevation_tifs[i])
    elevation = rast(all_elevation_tifs[i])

    max_depth_union = max_stage_union - elevation

    output_filename = paste0(output_dir, '/', gsub('elevation0', 'marine_warning_depth_maxima', basename(all_elevation_tifs[i])))

    #write_stars(max_depth_union, output_filename, options=c('COMPRESS=DEFLATE'), NA_value=-3.39999999999999996e+38)
    writeRaster(max_depth_union, output_filename, gdal=c('COMPRESS=DEFLATE'), NAflag=-3.39999999999999996e+38, overwrite=TRUE)

    return(TRUE)
}

# Run in parallel in this function, so failures will not cause all jobs to fail.
parallel_fun<-function(i){
    try(process_file_i(i))
}

# Run it.
library(parallel)
mclapply(1:length(all_marine_warning_tifs), parallel_fun, mc.cores=12, mc.preschedule=TRUE)
