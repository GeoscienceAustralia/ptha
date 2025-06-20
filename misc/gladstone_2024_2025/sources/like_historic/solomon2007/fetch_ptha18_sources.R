library(terra)
library(rptha)

cwd <- getwd()
setwd('/g/data/w85/tsunami/CODE/gadi/ptha_mm/ptha_access')
source('get_PTHA_results.R')


get_closest_sources<-function(name='solomon2', lon=156.34, lat=-7.79, n=4) {
    solomon = get_source_zone_events_data(name)
    solomon_sources = solomon$unit_source_statistics

    # Calculate the distance to the given point for each row
    solomon_sources$distance <- sqrt((solomon_sources$lon_c - lon)^2 + (solomon_sources$lat_c - lat)^2)

    # get the n closest unit sources
    closest_rows = list()
    for (i in 1:n) {
        # get the closest unit source
        closest_row <- solomon_sources[which.min(solomon_sources$distance), ]
        closest_rows[[i]] = closest_row
        # remove it from the solomon list
        solomon_sources <- solomon_sources[-which.min(solomon_sources$distance), ]
    }

    return(closest_rows)
}

get_tif_name <- function(unit_source) {
    return(basename(unit_source$initial_condition_file))
}

scale_unit_source <- function(unit_sources, Mw, dir_read, dir_write, shear_modulus=3.0e10, filter=TRUE) {
    for (source in unit_sources) {
        area = source$length * 1000 * source$width * 1000
        M0 = M0_2_Mw(Mw, inverse=TRUE)
        alpha = M0 / (shear_modulus * area)

        # Get the name of the tif file
        tif_name = paste0(dir_read, '/', get_tif_name(source))

        # Read the tif file
        tif = rast(tif_name)

        # Scale the tif file
        tif = tif * alpha

        path_write = paste0(dir_write, '/', gsub('.tif', '', get_tif_name(source)), '_alpha_', round(alpha, 3), '.tif')
        
        if (filter) {
            tif = kajiura_filter_source(tif, source$depth * 1000)
            path_write = gsub('.tif', '_filtered.tif', path_write)
        }

        writeRaster(tif, path_write, overwrite=TRUE)
    }
    return(unit_sources)
}

kajiura_filter_source <- function(tif, depth) {
    # project to local cartesian
    tif <- project(tif, '+proj=utm +zone=45 +south +ellps=WGS84')
    xyz_df <- as.data.frame(tif, xy=TRUE)
    names(xyz_df) <- c('x', 'y', 'z')
    xyz <- xyz_df[, c('x', 'y', 'z')]
    xyz_filtered <- kajiura_filter(xyz, depth)
    tif_filtered <- rast(xyz_filtered, tif)
    return(tif_filtered)
}


closest_sources = get_closest_sources('solomon2', 156.34, -7.79, n=4)
setwd(cwd)
warning('The unit sources were copied into the dir_read directory from the SOURCE_ZONES in the hazard directory. You will need to the same or script this up.')
dir_read = 'ptha_unit_sources'
dir_write = 'ptha_unit_sources_scaled'
dir.create(dir_write)
closest_sources <- scale_unit_source(closest_sources, Mw=8.1, dir_read=dir_read, dir_write=dir_write, filter=FALSE)
