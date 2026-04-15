#
# Make lines defining "breakwalls" (i.e. features that we'd like to ensure are
# burned into the elevation grid, irrespective of its resolution).
#
# The idea is the lines closely track the local elevation maxima representing
# barriers to the flow. Our hydrodynamic model grids are not always fine enough
# to be sure to represent these properly without this treatment.
#
# Given an input line shapefile showing where the breakwall is, we interpolate
# along the line with a prescribed point spacing. At each interpolated point,
# we search an associated elevation raster for the max-elevation within
# some buffer-radius. That max-elevation is associated with the interpolated
# point. Finally the resulting x/y/z points are exported to csv (and later read
# into our model and burned into the DEM). 
#
# The use of the 'buffer-radius' helps prevent us accidently using a lower
# elevation value (e.g. because the DEM cells do not perfectly track the
# elevation max, or because of the details regarding how these values are later
# burned into the model elevation). The values are set by judgment + inspection of
# the effect on the model.
#

library(rptha)
source('convert_lines_to_lon_lat_maxelev.R')

# Workhorse function here -- for all shapefiles associated with a dem
# (local_raster), make a xyz csv file containing the line (with a spacing of
# the provided distance in m). The z values contain the MAXIMUM cell in a
# buffer of size line_spacing around each point.
make_lines_with_max_buffer<-function(all_shp, local_raster, line_spacing, buffer_size, add_to_z=0){

    all_lines = vector(mode='list', length=length(all_shp))

    # Loop over each shapefile and convert to x/y/z csv file
    for(i in 1:length(all_shp)){
        all_lines[[i]] = convert_lines_to_lon_lat_maxelev(all_shp[i], local_raster, 
                                                          line_spacing_m = line_spacing,
                                                          buffer_width_m = buffer_size)
        colnames(all_lines[[i]][[1]]) = c('lon', 'lat', 'z')
    
        # At Carnarvon, the onshore Lidar has a significant offset relative to
        # the other data, so include that.
        all_lines[[i]][[1]][,3] = all_lines[[i]][[1]][,3] + add_to_z

        write.csv(all_lines[[i]][[1]], file=gsub('shp', 'csv', all_shp[i]), row.names=FALSE, 
                  quote=FALSE)
    }
}

make_breakwall <- function(shape_files, raster_file, line_spacing=3, buffer_size=5, add_to_z=0){
    shape_files = Sys.glob(shape_files)
    # Point elevation based on maxima inside this radius
    make_lines_with_max_buffer(shape_files, raster_file, line_spacing, buffer_size, add_to_z)
}

ALL_RASTERS = readLines('../../elevation/swals_elevation_files_in_preference_order.txt')

#
# Site-specific applications below here
#

# Geraldton
geraldton_horrocks_raster = ALL_RASTERS[grep("Area_B-C-D-E-F_Coastal.tif", ALL_RASTERS)]
stopifnot(length(geraldton_horrocks_raster) == 1)
make_breakwall("geraldton/*.shp", geraldton_horrocks_raster)

# Useless Loop
bathy_lidar_raster = ALL_RASTERS[grep("all_WADoT_2025_bathytopo_AHD_WGS84.vrt", ALL_RASTERS)]
stopifnot(length(bathy_lidar_raster) == 1)
make_breakwall("uselessloop/*.shp", bathy_lidar_raster)

# Carnarvon
# Use the onshore lidar BUT subtract 0.34m from it for consistency with the onshore Lidar.
# The offset is visually obvious when comparing the two datasets.
# Richard Nicholls confirmed both datasets with LANDGATE survey marks, finding the onshore
# lidar had the larger of the errors.
# See PDF of those email trail stored alongside the onshore lidar.
onshore_lidar_carnarvon = ALL_RASTERS[grep("all_WADoT_onshore_lidar_tifs_clipped.vrt", ALL_RASTERS)]
offset_carnarvon = -(0.079 - (-1*0.258))
print(paste0('WARNING: Applying offset of ', offset_carnarvon, ' to Carnarvon breakwalls using the onshore lidar ', onshore_lidar_carnarvon))
make_breakwall('carnarvon/*.shp', onshore_lidar_carnarvon, add_to_z=offset_carnarvon)

#
# Manual walls
#

## Interpolate between p1, p2
#interpolate_points <- function(p1, p2, num_points) {
#  wall <- matrix(NA, ncol = 3, nrow = num_points)
#  for (i in 1:num_points) {
#    wall[i, ] <- ((i - 1) * p2 + (num_points - i) * p1) / (num_points - 1)
#  }
#  return(wall)
#}
#
## Make manual floodgate
#make_manual_floodgate <- function(p1, p2, dir, filename) {
#    num_points = 5
#    wall = interpolate_points(p1, p2, num_points)
#    wall = as.data.frame(wall)
#    names(wall) = c('lon', 'lat', 'z')
#    dir.create(dir, showWarnings=FALSE)
#    pathwrite = paste0(dir, "/", "filename")
#    write.csv(wall, pathwrite, row.names=FALSE, quote=FALSE)
#}
