# Make some xyz lines defining breakwalls, with elevations set manually based on the final part of the directory name
library(rptha)

all_shp = Sys.glob('*/*.shp')
line_spacing_m = 3

# Ensure the output filenames will be unique
all_output_csv = gsub('.shp', '.csv', basename(all_shp))
if(length(unique(all_output_csv)) != length(all_output_csv)){
    stop('Shapefile names must be unique to avoid over-writing files')
}

for(i in 1:length(all_shp)){

    line_shapefile = all_shp[i]

    # Read lines denoting approximate levee positions
    lines = readOGR(line_shapefile, gsub('.shp', '', basename(line_shapefile)))

    if(!isLonLat(lines)) stop('line_shapefile must be in lon-lat coordinate system')

    # Interpolate with prescribed spacing
    lines_interp = approxSpatialLines(lines, spacing=line_spacing_m/1000, longlat=TRUE, 
                                      distinguish_disjoint_line_segments=TRUE)

    # Set elevation (from filename). Take the final part (after '_') and make numeric
    dirsplit = strsplit(dirname(line_shapefile), '_')[[1]]
    line_elev = as.numeric(dirsplit[length(dirsplit)])
    stopifnot(is.finite(line_elev))
    lines_xyz = cbind(coordinates(lines_interp), rep(line_elev, length(lines_interp)))

    colnames(lines_xyz) = c('lon', 'lat', 'z')

    # Write output
    output_file = gsub('.shp', '.csv', basename(line_shapefile), fixed=TRUE)

    write.csv(lines_xyz, file=output_file, row.names=FALSE)
}
