#
# Contour the arrival time rasters
#

library(terra)

arrival_time_dirs = Sys.glob('../NSW_tsunami_modelling_project_final_outputs/Arrival_times/*/*')
contour_levels = seq(0, 40*3600, by=300) # 300 second contour interval over 40 hours

# Loop over the individual directories (each with a single vrt)
for(i in 1:length(arrival_time_dirs)){

    arrival_time_dir_i = arrival_time_dirs[i]
    source_zone = basename(dirname(arrival_time_dir_i))
    arrival_type = basename(arrival_time_dir_i)

    print(paste0('Working on ', source_zone, ' ', arrival_type))

    # Make directories for the outputs -- store the contours of each individual tile separately
    output_dir = paste0(source_zone, '/', arrival_type)
    dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
    output_dir_individual = paste0(output_dir, '/individual_tiles')
    dir.create(output_dir_individual, recursive=TRUE, showWarnings=FALSE)

    # Find the tif filenames in the directory
    all_rasts = Sys.glob(paste0(arrival_time_dir_i, '/*.tif'))

    # Make a contour for each tif separately -- much faster than contouring the vrt.
    for(j in 1:length(all_rasts)){
        raster_j = all_rasts[j]
        rj = rast(raster_j)
        rj_contour = try(as.contour(rj, levels=contour_levels, maxcells=1e+07))
        # Contouring could fail if the raster as only missing-data pixels
        if(!is(rj_contour, 'try-error')){
            # It worked - save the file
            output_file = paste0(output_dir_individual, '/', gsub('.tif', '.shp', basename(raster_j), fixed=TRUE))
            writeVector(rj_contour, output_file, filetype='ESRI Shapefile', overwrite=TRUE)
        }
    }

    # Merge the individual contours into a single file
    tile_contour_files = Sys.glob(paste0(output_dir_individual, '/*.shp'))
    tile_contours = lapply(tile_contour_files, vect)
    merged_tile_contours = do.call(rbind, tile_contours)
    writeVector(merged_tile_contours, 
        paste0(output_dir, '/', source_zone, '_', arrival_type, '_merged_contours.shp'),
        filetype = 'ESRI Shapefile',
        overwrite=TRUE)
}
