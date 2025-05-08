# Given a line_shapefile, a dem, and a line_spacing in meters, interpolate along the line with
# the prescribed spacing and find the associated elevation value from the dem. APPLY A MAX-BUFFER
# AROUND EACH POINT WITH RADIUS = BUFFER_WIDTH_M
convert_lines_to_lon_lat_maxelev<-function(line_shapefile, dem_tif_file, line_spacing_m, buffer_width_m){
    library(rptha)
    library(parallel)

    # Read lines denoting approximate levee positions
    lines = readOGR(line_shapefile, gsub('.shp', '', basename(line_shapefile)))

    if(!isLonLat(lines)) stop('line_shapefile must be in lon-lat coordinate system')

    # Interpolate with prescribed spacing
    lines_interp = approxSpatialLines(lines, spacing=line_spacing_m/1000, longlat=TRUE, 
                                      distinguish_disjoint_line_segments=TRUE)

    # Get the maximum DEM elevation on each interpolated point, within some buffer distance
    dem = raster(dem_tif_file)
    ncores = 1 # Used to use parallel to to the extraction. Not now. 
    my_inds = splitIndices(length(lines_interp[,1]), ncores)
    lines_z_vals = mclapply(my_inds, 
                            f<-function(inds) extract(dem, lines_interp[inds,], buffer=buffer_width_m, fun=max),
                            mc.cores=ncores)
    lines_all_z = unlist(lines_z_vals)

    # Get unique identifier for each line
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
