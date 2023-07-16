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
local_raster = ALL_RASTERS[grep("tile_6.tif", ALL_RASTERS, fixed=TRUE)]
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

