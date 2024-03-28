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
# point. (For some other breakwalls, we set these elevations manually, see ./MANUAL_WALLS).
# Finally the resulting x/y/z points are exported to csv (and later read
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
all_elevation_files = readLines('../elevation/swals_elevation_files_in_preference_order.txt')

# NSW LADS DEM
bathy2018 = all_elevation_files[grep('nsw_bathy2018_patched.vrt', all_elevation_files)]

# Find all shapefiles one directory down from the current directory
all_shp = Sys.glob('*/*.shp')

if(any(grepl('MANUAL_WALLS', all_shp))){
    print('ERROR: The script is finding a shapefile inside the MANUAL_WALLS directory.')
    print('       That location is used for walls that are set by other means (not using this script).')
    print('       Either move shapefiles in MANUAL_WALLS to a subdirectory in this folder (if their elevation is set with the DEM,')
    print('           or move shapefiles in MANUAL_WALLS to a subdirectory of that folder (if their elevation is set separately).')
}

#
# SET WALL ELEVATIONS WITH THE DEM
#

local_raster = bathy2018
# Points with this spacing along the line
line_spacing = 3
# Point elevation based on maxima inside this radius
buffer_size = 15
make_lines_with_max_buffer(all_shp, local_raster, line_spacing, buffer_size)

#
# SET WALL ELEVATIONS MANUALLY
#
local_env = new.env()
source('MANUAL_WALLS/make_manual_walls.R', local=local_env, chdir=TRUE)

#
# Make file list
#
breakwalls_files = Sys.glob('*/*.csv')
writeLines(breakwalls_files, 'breakwalls_file_list_relative_path.txt')
