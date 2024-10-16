#' Make csv with fixed elevation
#' 
#' Make a csv file containing the lon, lat, z of points along the line.
#' Lines are specified in shapefiles below.
library(rptha)


convert_lines_to_lon_lat <- function(line_shapefile, line_spacing_m) {
    # Read lines denoting approximate levee positions
    lines = readOGR(line_shapefile, gsub('.shp', '', basename(line_shapefile)))
    if(!isLonLat(lines)) stop('line_shapefile must be in lon-lat coordinate system')
    # Interpolate with prescribed spacing
    lines_interp = approxSpatialLines(lines, spacing=5/1000, longlat=TRUE, distinguish_disjoint_line_segments=TRUE)
    # Get the coordinates of these points
    line_coords = coordinates(lines_interp)
    return(line_coords)
}

#' make a csv file with given elevation. 
#' 
#' Containing lat, lon and z, where z is a given scalar elevation value
make_csv_with_elev<-function(shape_file, z, line_spacing=5){
    line_coords = convert_lines_to_lon_lat(shape_file, line_spacing_m=line_spacing)
    colnames(line_coords) = c('lon', 'lat')
    # add z column with fixed values
    z_df <- data.frame(z = rep(z, times=nrow(line_coords)))
    xyz_df <- cbind(line_coords, z_df)
    write.csv(xyz_df, file=gsub('shp', 'csv', shape_file), row.names=FALSE, quote=FALSE)
}

### Generate a csv file for each shape file below
# Geraldton
make_csv_with_elev("shapes/geraldton_champion_bay_n.shp", z=2.3)
make_csv_with_elev("shapes/geraldton_champion_bay_s.shp", z=2.0)
make_csv_with_elev("shapes/geraldton_wharf.shp", z=2.6)
# Jurian bay
make_csv_with_elev("shapes/jurian_bay_1.shp", z=2.7+0.6)  # wall is 60cm high
make_csv_with_elev("shapes/jurian_bay_2.shp", z=2.8+0.6)  # wall is 60cm high
make_csv_with_elev("shapes/jurian_bay_3.shp", z=2.7+0.6)  # wall is 60cm high
