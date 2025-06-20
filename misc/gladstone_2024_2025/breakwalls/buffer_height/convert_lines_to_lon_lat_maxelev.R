# Given a line_shapefile, a dem, and a line_spacing in meters, interpolate along the line with
# the prescribed spacing and find the associated elevation value from the dem. By default apply a max-buffer
# around each point with radius = buffer_width_m. We can convert that to a min-buffer with the optional 
# aggregation_fun argument.
library(rptha)
# library(parallel)


convert_lines_to_lon_lat_maxelev<-function(
    line_shapefile,
    dem_tif_file,
    line_spacing_m,
    buffer_width_m,
    aggregation_fun=max){

    # Read lines denoting approximate levee positions
    lines = readOGR(line_shapefile, gsub('.shp', '', basename(line_shapefile)))

    if(!isLonLat(lines)) stop('line_shapefile must be in lon-lat coordinate system')

    # Interpolate with prescribed spacing
    print("Interpolating lines")
    lines_interp = approxSpatialLines(
        lines,
        spacing=line_spacing_m/1000,
        longlat=TRUE,
        distinguish_disjoint_line_segments=TRUE
        )

    # Get the max DEM elevation on each interpolated point, within some buffer distance (or
    # use another aggregation function to get some other value)
    print("Reading DEM - possibly large memory required if using a large raster.")
    dem = raster(dem_tif_file)
    lines_all_z = extract(dem, lines_interp, buffer = buffer_width_m, fun=aggregation_fun)

    # Get unique identifier for each line
    print("Getting unique lines")
    unique_lines = unique(lines_interp@data[,1])

    # Pack into a list of x/y/z coordinates for output
    ul = list()
    clines = coordinates(lines_interp)
    for(i in 1:length(unique_lines)){
        k = which(lines_interp@data[,1] == unique_lines[i])
        ul[[i]] = cbind(clines[k,1], clines[k,2], lines_all_z[k])
    }

    return(ul)
}
