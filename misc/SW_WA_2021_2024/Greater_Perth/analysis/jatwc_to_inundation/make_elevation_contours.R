#
# Make elevation contours in a given warning zone
#
source('make_vrt.R')
library(terra)

# Use the model PTHA elevation rasters -- note only a subset of these will be in our warning zone.
ptha_elevation_rasts = Sys.glob('../../swals/OUTPUTS/Fuji_andaman2004_24hrs_domain010322_timevaryingRealistic-full-ambient_sea_level_0.0/RUN_20220303_092730617/elevation0*.tif')

# Get a set of rasters covering the zone -- here we choose no_threat, although any choice would be fine.
# The purpose is just to easily identify which domains should be included.
ATWS_ZONE_NAME = commandArgs(trailingOnly=TRUE)[1]
ATWS_Zone_name_nospace = gsub(' ', '-', ATWS_ZONE_NAME)
zone_files_no_threat = Sys.glob(paste0('Inundation_zones/', ATWS_Zone_name_nospace, '/no_threat*.tif'))

# Get files in the PTHA raster that match these files. 
# Not all are included, because the JATWC zone will not use all the model tifs
zone_files_endofname = sapply(zone_files_no_threat, function(x) {y = strsplit(x, '_')[[1]]; return(y[length(y)])})
ptha_rasts_endofname = sapply(ptha_elevation_rasts, function(x) {y = strsplit(x, '_')[[1]]; return(y[length(y)])})
matching_files = match(ptha_rasts_endofname, zone_files_endofname)

# Check we got a 1:1 match
stopifnot(all(sort(na.omit(matching_files)) == (1:length(zone_files_no_threat))  ))
k = which(!is.na(matching_files))
ptha_elevation_vrt = rast(make_vrt(ptha_elevation_rasts[k]))

# If the following parameter is too small, then some contours may be missing (and we get warnings).
# Must be less than .Machine$integer.max
options("max.contour.segments" = 2e+09) 
elevation_contours = as.contour(ptha_elevation_vrt, maxcells=Inf, levels=seq(0, 10, by=1))

# Write to shapefile
output_dir = paste0('elevation_contours/', ATWS_Zone_name_nospace, '_contours')
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
writeVector(elevation_contours, file=paste0(output_dir, '/', basename(output_dir), '.shp'), overwrite=TRUE)

#
# ALTERNATIVE
# library(stars)
# ptha_elevation_vrt = read_stars(st_mosaic(ptha_elevation_rasts[k], options=c('-resolution', 'highest')))
# elevation_contours = st_contour(ptha_elevation_vrt)
# st_write(elevation_contours, dsn='test_contours', layer='test_contours.shp', driver='ESRI Shapefile')
#
