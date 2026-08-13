#
# Make rasters where each pixel has the smallest JATWC-H value in each coastal
# zone that causes flooding. Run with (e.g.):
#
#     Rscript compute_scenario_statistics_in_zone.R 'Bunbury Geographe Coast'
#     Rscript map_threat_levels_in_zone.R 'Bunbury Geographe Coast'
#

library(sf)
library(units)
library(raster)

#
# INPUTS
#
asfm = new.env()
source('application_specific_metadata.R', local=asfm)

ATWS_ZONE_NAME = commandArgs(trailingOnly=TRUE)[1] # 'Perth Coast'

FLOODING_DEPTH = 0.001 # Depths above this are treated as flooded
print(paste0('Treating cells as flooded if the depth is >= ', FLOODING_DEPTH))

#
# END INPUTS
#

# ATWS Zones
atws_zones = asfm$atws_zones 
warning_zone = atws_zones[which(atws_zones$ATWS_Zones == ATWS_ZONE_NAME),]
if(nrow(warning_zone) != 1) stop('ATWS_ZONE_NAME should match exactly one value in attribute ATWS_Zones')

is_offshore_island = any(ATWS_ZONE_NAME == asfm$JATWC_offshore_island_zones)

# How many cores to use?
MC_CORES = asfm$DEFAULT_MC_CORES

# Determine which rasters from the multidomain are touching the warning zone
get_max_depth_rasters_touching_warning_zone<-function(
    scenario_raster_tar, 
    warning_zone, 
    STARTING_DIR,
    warning_zone_buffer_degrees){

    # Get filenames without extraction for speed
    all_raster_basenames = system(paste0('tar --list -f ', scenario_raster_tar), intern=TRUE)
    all_tifs = paste0('/vsitar/', scenario_raster_tar, '/',
        all_raster_basenames[grep('depth_as_max_stage_minus_elevation0_', all_raster_basenames)])

    overlapping = rep(NA, length(all_tifs))

    # Generous buffer in warning zone to include tifs
    using_s2 = sf_use_s2() # workaround for problems with sf::st_buffer, which makes jaggedy edges in buffering
    sf_use_s2(FALSE) # Force use of GEOS, which treats coordinates like cartesian coordinates
    warning_zone_buffer = st_buffer(warning_zone, 
        dist=set_units(warning_zone_buffer_degrees, 'arc_degree') 
            # Above buffer dist in degrees. Because we are using GEOS the units don't seem to matter.
        )
    sf_use_s2(using_s2)

    # For each raster, figure out if the extent touches the warning zone
    #for(i in 1:length(all_tifs))
    parfun<-function(i){
        r1_tiff = raster(all_tifs[i])
        r1_dx = res(r1_tiff)

        if(any(r1_dx > asfm$skip_if_cellsize_above_threshold) ){
            overlapping = FALSE
        }else{
            r1 = extent(r1_tiff)
            r1_bbox = matrix(
                c(r1@xmin, r1@ymin, 
                  r1@xmax, r1@ymin, 
                  r1@xmax, r1@ymax, 
                  r1@xmin, r1@ymax, 
                  r1@xmin, r1@ymin), ncol=2, byrow=TRUE)
            overlapping = (st_intersects(warning_zone_buffer, st_polygon(list(r1_bbox)), sparse=FALSE)[1,1])
        }
        return(overlapping)
    }
    library(parallel)
    is_overlapping = mclapply(1:length(all_tifs), parfun, mc.cores=MC_CORES, mc.preschedule=FALSE)
    overlapping = which(unlist(is_overlapping))

    overlapping_tifs = basename(all_tifs[overlapping])

    return(overlapping_tifs)
}

# Given a raster tile (single domain in the SWALS model), make a raster
# containing the minimum JATWC_H statistic causing flooding at each pixel.
get_min_H_causing_flooding<-function(
    raster_tile_name, 
    all_scenario_JATWC_H){

    # For each scenarios, domain rasters are stored in tar archives -- GDAL can
    # read inside tar archives using the right file name notation ('/vsitar/')
    all_raster_files = paste0('/vsitar/', 
        unlist(lapply(all_scenario_JATWC_H, function(y) y$model_dir)), 
        '/raster_output_files.tar/', raster_tile_name)

    # The rasters use missing data to denote dry areas. To prevent these wrecking the minima
    # computation, replace NA with a large positive number (larger than the largest JATWC_H value AND than the largest raster value)
    largest_H = max(unlist(lapply(all_scenario_JATWC_H, function(y) y$jatwc_H)))
    MISSING_VALUE =  9999999
    if(largest_H >= MISSING_VALUE){
        print(paste0('Unexpectedly large JATWC_H value: ', largest_H))
        MISSING_VALUE = MISSING_VALUE + largest_H
    }

    raster_extrema = raster(all_raster_files[1])
    raster_extrema = setValues(raster_extrema, MISSING_VALUE)

    raster_extrema_mat = as.matrix(raster_extrema)

    # Maxima computation in a loop
    for(i in 1:length(all_raster_files)){

        # Maxima computation is faster using matrices, rather than rasters
        # Consider stars -- I think like `t(read_stars(all_raster_files[i])[[1]])`
        r1 = as.matrix(raster(all_raster_files[i]))
        stopifnot(all(dim(r1) == dim(raster_extrema_mat)))

        # Fix for imperfectly recognized missing data values
        r1[r1 < asfm$raster_na_below] = NA
        r1[is.na(r1)] = MISSING_VALUE

        # Check we don't have pixels larger than MISSING VALUE (which would break the logic)
        max_r1 = max(r1)
        if( (!is.finite(max_r1)) | (max_r1 > MISSING_VALUE)){
            stop(paste0('Unexpectedly large depth value (breaking the code logic): ', 
                max_r1, ' ', all_raster_files[i]))
        }

        # We do not want cells with depth below FLOODING_DEPTH to be treated as flooded,
        # so instead treat them as missing
        r1[r1 < FLOODING_DEPTH] = MISSING_VALUE

        # All remaining non-missing pixels had depth values exceeding
        # FLOODING_DEPTH. Set these equal to the JATWC H value
        r1[(r1 < MISSING_VALUE)] = all_scenario_JATWC_H[[i]]$jatwc_H

        # Minimum JATWC_H
        raster_extrema_mat = pmin(raster_extrema_mat, r1)

    }

    rm(r1); gc()

    # Consider stars::st_as_stars
    raster_extrema = setValues(raster_extrema, raster_extrema_mat)

    # Anything that remains missing should really be NA
    raster_extrema[raster_extrema == MISSING_VALUE] = NA

    return(raster_extrema)
}

#
# Main program here
#

ATWS_Zone_name_nospace = gsub(' ', '-', ATWS_ZONE_NAME)

## Make space for outputs
working_dir = paste0('Inundation_zones/', ATWS_Zone_name_nospace)
dir.create(working_dir, showWarnings=FALSE, recursive=TRUE)
setwd(working_dir)

# Useful to return here
STARTING_DIR = getwd()

# Read values of JATWC_H that were computed with an earlier script
all_JATWC_H_file = Sys.glob('all-JATWC-H_*.RDS'); stopifnot(length(all_JATWC_H_file) == 1)

# Check the H_file name is consistent with our Zone name
stopifnot(grepl(ATWS_Zone_name_nospace, all_JATWC_H_file))
all_scenario_JATWC_H = readRDS(all_JATWC_H_file)

# Find all rasters that touch this warning zone (with some buffering of the
# warning zone so we get onshore areas)
max_depth_rasters_touching_warning_zone = get_max_depth_rasters_touching_warning_zone(
    asfm$all_scenario_raster_tars[1], warning_zone, STARTING_DIR, asfm$warning_zone_buffer_degrees)

# Make rasters for a given JATWC_H_ranges entry, on a single tile
tile_function<-function(raster_tile){
    library(sf)
    library(units)
    library(raster)

    min_jatwc_h_flooding = get_min_H_causing_flooding(raster_tile, all_scenario_JATWC_H)

    output_raster = paste0('min_jatwc_h_causing_flooding_', raster_tile)
    writeRaster(min_jatwc_h_flooding, output_raster, options=c('COMPRESS=DEFLATE'), overwrite=TRUE)

    rm(min_jatwc_h_flooding, output_raster); gc()
    return(0)
}


# Allow failures in parallel
parallel_fun_combined<-function(raster_tile) try(tile_function(raster_tile))


# Make all the rasters (using max-stage)
library(parallel)
my_cluster = makeCluster(MC_CORES)
export_job = clusterExport(my_cluster, varlist=ls())
# Combined max-stage and arrival time into one batch of work
combined_files = as.list(max_depth_rasters_touching_warning_zone)
all_results = parLapplyLB(my_cluster, combined_files, parallel_fun_combined, chunk.size=1)
any_try_errors = any(unlist(lapply(all_results, function(x) is(x, 'try-error'))))
stopCluster(my_cluster)

if(any_try_errors){
    print(all_results)
    stop('There were failures computing min_jatwc_h_causing_flooding')
}

