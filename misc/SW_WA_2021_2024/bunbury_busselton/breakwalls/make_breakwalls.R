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


#
# Site-specific applications below here
#
ALL_RASTERS = readLines('../elevation/swals_elevation_files_in_preference_order.txt')

# WA DoT LAS/LIDAR
file_home = '/media/gareth/Windows7_OS/Users/gareth/Documents/work/AustPTHA/DATA/ELEV/WA/WA_DoT_DSM/GreaterPerth_5m_WGS84.vrt'
file_NCI = '/g/data/w85/tsunami/DATA/ELEVATION/WA/WA_DoT_DSM/GreaterPerth_5m_WGS84.vrt'
wa_dot_las_lidar = ifelse(file.exists(file_NCI), file_NCI, file_home)
print('Using WA_DoT_LAS_LIDAR for walls near Fremantle boat harbour (sensitivity test)')


# Geraldton
all_shp = Sys.glob('perth/geraldto*.shp')
local_raster = ALL_RASTERS[grep("tile_32.tif", ALL_RASTERS, fixed=TRUE)]
line_spacing = 5
# Point elevation based on maxima inside this radius
buffer_size = 15
make_lines_with_max_buffer(all_shp, local_raster, line_spacing, buffer_size)

# JurianBay
all_shp = Sys.glob('perth/jurian*.shp')
local_raster = ALL_RASTERS[grep("tile_24.tif", ALL_RASTERS, fixed=TRUE)]
line_spacing = 5
# Point elevation based on maxima inside this radius
buffer_size = 15
make_lines_with_max_buffer(all_shp, local_raster, line_spacing, buffer_size)

# TwoRocks
all_shp = Sys.glob('perth/twoRock*.shp')
local_raster = ALL_RASTERS[grep("tile_1.tif", ALL_RASTERS, fixed=TRUE)]
line_spacing = 5
# Point elevation based on maxima inside this radius
buffer_size = 15
make_lines_with_max_buffer(all_shp, local_raster, line_spacing, buffer_size)

# North of Hillarys
all_shp = Sys.glob('perth/northP*.shp')
local_raster = ALL_RASTERS[grep("tile_2.tif", ALL_RASTERS, fixed=TRUE)]
line_spacing = 5
# Point elevation based on maxima inside this radius
buffer_size = 15
make_lines_with_max_buffer(all_shp, local_raster, line_spacing, buffer_size)


# Hillarys
all_shp = Sys.glob('perth/hillar*.shp')
local_raster = ALL_RASTERS[grep("tile_3.tif", ALL_RASTERS, fixed=TRUE)]
line_spacing = 5
# Point elevation based on maxima inside this radius
buffer_size = 15
make_lines_with_max_buffer(all_shp, local_raster, line_spacing, buffer_size)

# Near mouth of swan-canning estuary
all_shp = Sys.glob('perth/swanCannin*.shp')
#local_raster = ALL_RASTERS[grep("tile_6.tif", ALL_RASTERS, fixed=TRUE)]
local_raster = wa_dot_las_lidar
line_spacing = 5
# Point elevation based on maxima inside this radius
buffer_size = 15
make_lines_with_max_buffer(all_shp, local_raster, line_spacing, buffer_size)

# South of Perth, including on Garden Island
all_shp = Sys.glob('perth/southpert*.shp')
local_raster = ALL_RASTERS[grep("tile_7.tif", ALL_RASTERS, fixed=TRUE)]
line_spacing = 5
# Point elevation based on maxima inside this radius
buffer_size = 15
make_lines_with_max_buffer(all_shp, local_raster, line_spacing, buffer_size)

# Near Mandurah
all_shp = Sys.glob('perth/mandurah*.shp')
local_raster = ALL_RASTERS[grep("tile_9.tif", ALL_RASTERS, fixed=TRUE)]
line_spacing = 5
# Point elevation based on maxima inside this radius
buffer_size = 15
make_lines_with_max_buffer(all_shp, local_raster, line_spacing, buffer_size)

# Peppermint Grove
all_shp = Sys.glob('bunbury_busselton/peppermintGrove*.shp')
local_raster = ALL_RASTERS[grep("Mosaic1mDEM163_WGS84", ALL_RASTERS, fixed=TRUE)]
line_spacing = 5
# Point elevation based on maxima inside this radius
buffer_size = 10
make_lines_with_max_buffer(all_shp, local_raster, line_spacing, buffer_size)

# Bunbury
all_shp = Sys.glob('bunbury_busselton/bunburyWall*.shp')
local_raster = wa_dot_las_lidar
line_spacing = 5
# Point elevation based on maxima inside this radius
buffer_size = 15
make_lines_with_max_buffer(all_shp, local_raster, line_spacing, buffer_size)


# Other places around Busselton - Dunsborough
all_shp = Sys.glob('bunbury_busselton/*Buss*.shp')
#local_raster = ALL_RASTERS[grep("Mosaic1mDEM163_WGS84", ALL_RASTERS, fixed=TRUE)] # Shane Martin's DEM from the storm-surge study
local_raster = ALL_RASTERS[grep("Busselton_merged_20220330_Aerometrix_Lidar_1m_WGS84", ALL_RASTERS, fixed=TRUE)] # Newer DEM

line_spacing = 2
# Point elevation based on maxima inside this radius
buffer_size = 2
make_lines_with_max_buffer(all_shp, local_raster, line_spacing, buffer_size)


# A 2024 addition to Port Geographe that wasn't include in the original models
all_shp = Sys.glob('PortGeographUpdate/*.shp')
local_raster = ALL_RASTERS[grep("Busselton_merged_20220330_Aerometrix_Lidar_1m_WGS84", ALL_RASTERS, fixed=TRUE)] # Newer DEM
line_spacing = 2
# Point elevation based on maxima inside this radius
buffer_size = 2
make_lines_with_max_buffer(all_shp, local_raster, line_spacing, buffer_size)



#
# Manual flood-gates
#

print('Making 2 flood-gates at Wonnerup (manually)')

# Gate 1 -- elevation of 2m inferred from Shane Martin's high-res DEM used for 2014 storm surge study
p1 = rev(c(2, -33.61969776, 115.41416817))
p2 = rev(c(2, -33.619386  , 115.413722  ))
# Interpolate between p1, p2
wall = matrix(NA, ncol=3, nrow=5)
for(i in 1:5){
    wall[i,] = ( (i-1)*p2 + (5-i)*p1 )/4
}
wall = as.data.frame(wall)
names(wall) = c('lon', 'lat', 'z')
dir.create('wonnerup_floodgate', showWarnings=FALSE)
write.csv(wall, 'wonnerup_floodgate/wonnerup_floodgate_1.csv', row.names=FALSE, quote=FALSE)

# Gate 2 -- elevation of 2 m inferred from nearby topography
p1 = rev(c(2, -33.613692878, 115.425668470))
p2 = rev(c(2, -33.61331    , 115.426001   ))
# Interpolate between p1, p2
wall = matrix(NA, ncol=3, nrow=5)
for(i in 1:5){
    wall[i,] = ( (i-1)*p2 + (5-i)*p1 )/4
}
wall = as.data.frame(wall)
names(wall) = c('lon', 'lat', 'z')
dir.create('wonnerup_floodgate', showWarnings=FALSE)
write.csv(wall, 'wonnerup_floodgate/wonnerup_floodgate_2.csv', row.names=FALSE, quote=FALSE)



print('Making Bunbury floodgate "the plug" (manually)')
# Bunbury floodgate
p1 = c(115.64141925, -33.32074460, 2.26)
p2 = c(115.64153596, -33.32075600, 2.26)
# Interpolate between p1, p2
wall = matrix(NA, ncol=3, nrow=5)
for(i in 1:5){
    wall[i,] = ( (i-1)*p2 + (5-i)*p1 )/4
}
wall = as.data.frame(wall)
names(wall) = c('lon', 'lat', 'z')
dir.create('bunbury_floodgate', showWarnings=FALSE)
write.csv(wall, 'bunbury_floodgate/bunbury_floodgate.csv', row.names=FALSE, quote=FALSE)


