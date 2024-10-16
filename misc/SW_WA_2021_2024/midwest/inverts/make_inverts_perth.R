#
# Make lines defining "inverts" (i.e. features that we'd like to ensure are
# burned into the elevation grid as a min elevation, irrespective of its resolution).
#

library(rptha)
source('../breakwalls/convert_lines_to_lon_lat_maxelev.R')

# Workhorse function here -- for all shapefiles associated with a dem
# (local_raster), make a xyz csv file containing the line (with a spacing of
# the provided distance in m). The z values contain the MINIMUM cell in a
# buffer of size line_spacing around each point.
make_lines_with_min_buffer<-function(all_shp, local_raster, line_spacing, buffer_size){

    all_lines = vector(mode='list', length=length(all_shp))

    # Loop over each shapefile and convert to x/y/z csv file
    for(i in 1:length(all_shp)){
        all_lines[[i]] = convert_lines_to_lon_lat_maxelev(all_shp[i], local_raster, 
                                                          line_spacing_m = line_spacing,
                                                          buffer_width_m = buffer_size,
                                                          aggregation_fun = min)
        colnames(all_lines[[i]][[1]]) = c('lon', 'lat', 'z')
        write.csv(all_lines[[i]][[1]], file=gsub('shp', 'csv', all_shp[i]), row.names=FALSE, 
                  quote=FALSE)
    }
}

stop('The example below was from NSW -- currently not using inverts in this model')
#
# Site-specific applications below here
#
all_elevation_files = readLines('../elevation/swals_elevation_files_in_preference_order.txt')

# NSW LADS DEM
bathy2018 = all_elevation_files[grep('nsw_bathy2018_patched.vrt', all_elevation_files)]

# Shapefiles
all_shp = Sys.glob('*/*.shp')

local_raster = bathy2018
# Points with this spacing along the line
line_spacing = 3
# Point elevation based on minima inside this radius
buffer_size = 15
make_lines_with_min_buffer(all_shp, local_raster, line_spacing, buffer_size)

inverts_files = Sys.glob('*/*.csv')
writeLines(inverts_files, 'inverts_file_list_relative_path.txt')
