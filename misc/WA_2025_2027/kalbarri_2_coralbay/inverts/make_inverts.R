#
# Make lines defining "inverts" (i.e. features that we'd like to ensure are
# burned into the elevation grid as a min elevation, irrespective of its resolution).
#

library(rptha)
source('convert_lines_to_lon_lat_maxelev.R')

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
        write.csv(all_lines[[i]][[1]], file=gsub('shp', 'csv', all_shp[i]), row.names=FALSE, quote=FALSE)
    }
}

convert_lines_to_lon_lat <- function(line_shapefile, line_spacing_m) {
    # Read lines denoting approximate levee positions
    lines = readOGR(line_shapefile, gsub('.shp', '', basename(line_shapefile)))
    if(!isLonLat(lines)) stop('line_shapefile must be in lon-lat coordinate system')
    # Interpolate with prescribed spacing
    lines_interp = approxSpatialLines(lines, spacing=line_spacing_m/1000, longlat=TRUE, distinguish_disjoint_line_segments=TRUE)
    # Get the coordinates of these points
    line_coords = coordinates(lines_interp)
    return(line_coords)
}

#' make a csv file with given elevation. 
#' 
#' Containing lat, lon and z, where z is a given scalar elevation value
make_csv_with_elev<-function(shape_file, z, line_spacing=3){
    line_coords = convert_lines_to_lon_lat(shape_file, line_spacing_m=line_spacing)
    colnames(line_coords) = c('lon', 'lat')
    # add z column with fixed values
    z_df <- data.frame(z = rep(z, times=nrow(line_coords)))
    xyz_df <- cbind(line_coords, z_df)
    write.csv(xyz_df, file=gsub('shp', 'csv', shape_file), row.names=FALSE, quote=FALSE)
}

# Make fixed height inverts
make_csv_with_elev("fixed_height_lines/carnarvon1.shp", z=0.15)
make_csv_with_elev("fixed_height_lines/carnarvon2.shp", z=0.80)

# Add the files to a list for use in swals
invert_files <- Sys.glob('*/*.csv')
invert_files_full <- normalizePath(invert_files)
writeLines(invert_files_full, 'swals_invert_files.txt')
