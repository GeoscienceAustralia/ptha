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
make_lines_with_max_buffer<-function(all_shp, local_raster, line_spacing, buffer_size){

    all_lines = vector(mode='list', length=length(all_shp))

    # Loop over each shapefile and convert to x/y/z csv file
    for(i in 1:length(all_shp)){
        all_lines[[i]] = convert_lines_to_lon_lat_maxelev(all_shp[i], local_raster, 
                                                          line_spacing_m = line_spacing,
                                                          buffer_width_m = buffer_size)
        colnames(all_lines[[i]][[1]]) = c('lon', 'lat', 'z')
        write.csv(all_lines[[i]][[1]], file=gsub('shp', 'csv', all_shp[i]), row.names=FALSE, 
                  quote=FALSE)
    }
}

make_breakwall <- function(shape_files, raster_file, all_rasters, line_spacing=5, buffer_size=15) {
    shape_files = Sys.glob(shape_files)
    local_raster = all_rasters[grep(raster_file, all_rasters, fixed=TRUE)]
    # Point elevation based on maxima inside this radius
    make_lines_with_max_buffer(shape_files, local_raster, line_spacing, buffer_size)
}

#
# Site-specific applications below here
#
ALL_RASTERS = readLines('../../elevation/swals_elevation_files_in_preference_order.txt')

# WA DoT LAS/LIDAR
file_home = '/media/gareth/Windows7_OS/Users/gareth/Documents/work/AustPTHA/DATA/ELEV/WA/WA_DoT_DSM/GreaterPerth_5m_WGS84.vrt'
file_NCI = '/g/data/w85/tsunami/DATA/ELEVATION/WA/WA_DoT_DSM/GreaterPerth_5m_WGS84.vrt'
wa_dot_las_lidar = ifelse(file.exists(file_NCI), file_NCI, file_home)
print('Using WA_DoT_LAS_LIDAR for walls near Fremantle boat harbour (better captures breakwalls)')


# Midwest region modelling
make_breakwall("midwest/port_denison*.shp", "tile_29.tif", ALL_RASTERS)

# Perth region modelling
make_breakwall("perth/geraldto*.shp", "tile_32.tif", ALL_RASTERS)
make_breakwall("perth/jurian*.shp", "tile_24.tif", ALL_RASTERS)
make_breakwall("perth/twoRock*.shp", "tile_1.tif", ALL_RASTERS)
make_breakwall("perth/northP*.shp", "tile_2.tif", ALL_RASTERS)
make_breakwall("perth/hillar*.shp", "tile_3.tif", ALL_RASTERS)
# make_breakwall("perth/swanCannin*.shp", wa_dot_las_lidar, ALL_RASTERS)  # Error in swanCanning1?
make_breakwall("perth/southpert*.shp", "tile_7.tif", ALL_RASTERS)
make_breakwall("perth/mandurah*.shp", "tile_9.tif", ALL_RASTERS)

# Bunbury and Busselton region
make_breakwall("bunbury_busselton/peppermintGrove*.shp", "Mosaic1mDEM163_WGS84", ALL_RASTERS, buffer_size=10)
# note that wa_dot_las_lidar is not in ALL_RASTERS, so pass it as the rasters to grep through
make_breakwall("bunbury_busselton/bunburyWall*.shp", wa_dot_las_lidar, wa_dot_las_lidar)
# Other places around Busselton - Dunsborough
make_breakwall(
    "bunbury_busselton/*Buss*.shp",
    "Busselton_merged_20220330_Aerometrix_Lidar_1m_WGS84",
    ALL_RASTERS,
    line_spacing=2,
    buffer_size=2
    )


#
# Manual flood-gates
#

# Interpolate between p1, p2
interpolate_points <- function(p1, p2, num_points) {
  wall <- matrix(NA, ncol = 3, nrow = num_points)
  for (i in 1:num_points) {
    wall[i, ] <- ((i - 1) * p2 + (num_points - i) * p1) / (num_points - 1)
  }
  return(wall)
}

# Make manual floodgate
make_manual_floodgate <- function(p1, p2, dir, filename) {
    num_points = 5
    wall = interpolate_points(p1, p2, num_points)
    wall = as.data.frame(wall)
    names(wall) = c('lon', 'lat', 'z')
    dir.create(dir, showWarnings=FALSE)
    pathwrite = paste0(dir, "/", "filename")
    write.csv(wall, pathwrite, row.names=FALSE, quote=FALSE)
}

print('Making 2 flood-gates at Wonnerup (manually)')
# Gate 1 -- elevation of 2m inferred from Shane Martin's high-res DEM used for 2014 storm surge study
p1 = rev(c(2, -33.61969776, 115.41416817))
p2 = rev(c(2, -33.619386  , 115.413722  ))
make_manual_floodgate(p1, p2, "wonnerup_floodgate", "wonnerup_floodgate_1.csv")

# Gate 2 -- elevation of 2 m inferred from nearby topography
p1 = rev(c(2, -33.613692878, 115.425668470))
p2 = rev(c(2, -33.61331    , 115.426001   ))
make_manual_floodgate(p1, p2, "wonnerup_floodgate", "wonnerup_floodgate_2.csv")

# Bunbury floodgate
print('Making Bunbury floodgate "the plug" (manually)')
p1 = c(115.64141925, -33.32074460, 2.26)
p2 = c(115.64153596, -33.32075600, 2.26)
make_manual_floodgate(p1, p2, "bunbury_floodgate", "bunbury_floodgate.csv")
