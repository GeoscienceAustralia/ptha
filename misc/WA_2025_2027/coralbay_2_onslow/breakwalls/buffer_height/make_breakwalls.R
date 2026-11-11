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

# Exmouth
onshore_lidar_path <- "/g/data/w85/tsunami/DATA/ELEVATION/WA/WADoT_2025_new_lidar_onshore_only/tifs_wgs84_clipped_to_land/all_WADoT_onshore_lidar_tifs_clipped.vrt"
make_breakwall("shapes/exmouth_breakwater_1.shp", onshore_lidar_path)
make_breakwall("shapes/exmouth_breakwater_2.shp", onshore_lidar_path)
make_breakwall("shapes/exmouth_breakwater_3.shp", onshore_lidar_path)
make_breakwall("shapes/exmouth_breakwater_4.shp", onshore_lidar_path)

# Onslow salt sea walls
green_lidar_path <- "/g/data/w85/tsunami/DATA/ELEVATION/WA/WADoT_2025_new_lidar/all_WADoT_2025_bathytopo_AHD_WGS84_clean.vrt"
make_breakwall("shapes/onslow_salt_wall_1.shp", green_lidar_path)
