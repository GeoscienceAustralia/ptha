#
# Make rasters that correspond to JATWC zones. Run after the other script, like:
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

ATWS_ZONE_NAME = commandArgs(trailingOnly=TRUE)[1] # 'Capricornia Coast'

#
# END INPUTS
#

# ATWS Zones
atws_zones = asfm$atws_zones 
warning_zone = atws_zones[which(atws_zones$ATWS_Zones == ATWS_ZONE_NAME),]
if(nrow(warning_zone) != 1) stop('ATWS_ZONE_NAME should match exactly one value in attribute ATWS_Zones')

# How many cores to use?
MC_CORES = asfm$DEFAULT_MC_CORES

# Determine which rasters from the multidomain are touching the warning zone
get_max_stage_rasters_touching_warning_zone<-function(
    scenario_raster_tar, 
    warning_zone, 
    STARTING_DIR,
    warning_zone_buffer_degrees=asfm$warning_zone_buffer_degrees){

    # Get filenames without extraction for speed
    all_raster_basenames = system(paste0('tar --list -f ', scenario_raster_tar), intern=TRUE)
    all_tifs = paste0('/vsitar/', scenario_raster_tar, '/',
        all_raster_basenames[grep('max_stage_domain_', all_raster_basenames)])

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

# Compute the maxima (or minima) over all scenarios with JATWC parameter within a specified
# range, for one particular raster tile (i.e. one domain in the SWALS model).
get_raster_extremes_stratified_by_H<-function(
    raster_tile_name, 
    all_scenario_JATWC_H, 
    JATWC_H_range,
    operation='maxima'){

    # For each scenarios, domain rasters are stored in tar archives -- GDAL can
    # read inside tar archives using the right file name notation ('/vsitar/')
    all_raster_files = paste0('/vsitar/', 
        unlist(lapply(all_scenario_JATWC_H, function(y) y$model_dir)), 
        '/raster_output_files.tar/', raster_tile_name)

    if(operation == 'maxima'){
        # The rasters use missing data to denote dry areas. To prevent these wrecking the maxima
        # computation, replace NA with a large negative number
        MISSING_VALUE = -9999999
    }else if(operation == 'minima'){
        # The rasters use missing data to denote dry areas. To prevent these wrecking the minima
        # computation, replace NA with a large positive number
        MISSING_VALUE =  9999999
    }else{
        stop(paste0('unknown operation ', operation))
    }

    raster_extrema = raster(all_raster_files[1])
    raster_extrema = setValues(raster_extrema, MISSING_VALUE)

    raster_extrema_mat = as.matrix(raster_extrema)

    # Maxima computation in a loop
    for(i in 1:length(all_raster_files)){
        print(i)
        # Skip rasters that do not meet the H criteria
        if( (all_scenario_JATWC_H[[i]]$jatwc_H < JATWC_H_range[1]) | 
            (all_scenario_JATWC_H[[i]]$jatwc_H > JATWC_H_range[2]) ) next

        # Maxima computation is faster using matrices, rather than rasters
        # Consider stars -- I think like `t(read_stars(all_raster_files[i])[[1]])`
        r1 = as.matrix(raster(all_raster_files[i]))
        stopifnot(all(dim(r1) == dim(raster_extrema_mat)))

        # Fix for imperfectly recognized missing data values
        r1[r1 < asfm$raster_na_below] = NA
        r1[is.na(r1)] = MISSING_VALUE

        if(operation == 'maxima'){
            raster_extrema_mat = pmax(raster_extrema_mat, r1)
        }else if(operation == 'minima'){
            raster_extrema_mat = pmin(raster_extrema_mat, r1)
        }else{
            stop(paste0('unknown operation ', operation))
        }

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
max_stage_rasters_touching_warning_zone = get_max_stage_rasters_touching_warning_zone(
    asfm$all_scenario_raster_tars[1], warning_zone, STARTING_DIR)

# Make rasters for all JATWC_H_ranges on a single tile
tile_function<-function(raster_tile, operation, jatwc_h_ranges_index){
    library(sf)
    library(units)
    library(raster)

    #for(i in 1:length(asfm$JATWC_H_ranges)){
    i = jatwc_h_ranges_index

    domain_warning = get_raster_extremes_stratified_by_H(
        raster_tile, all_scenario_JATWC_H, asfm$JATWC_H_ranges[[i]], operation=operation)
    output_raster = paste0(names(asfm$JATWC_H_ranges)[i], '_', raster_tile)
    writeRaster(domain_warning, output_raster, options=c('COMPRESS=DEFLATE'), overwrite=TRUE)
    #}

    rm(domain_warning, output_raster); gc()
    return(0)
}


# This works for rasters with either max_stage OR arrival time, and
# lets us increase the parallel work
parallel_fun_combined<-function(combined_data_i){
    stopifnot(nrow(combined_data_i) == 1)
    stopifnot(all(names(combined_data_i) == c('raster_tile', 'jatwc_h_ranges_index')))
    
    raster_tile = combined_data_i$raster_tile
    jatwc_h_ranges_index = combined_data_i$jatwc_h_ranges_index

    if(grepl('max_stage', raster_tile)){
        x = try(tile_function(raster_tile, 'maxima', jatwc_h_ranges_index))
    }else if(grepl('arrival_time', raster_tile)){
        x = try(tile_function(raster_tile, 'minima', jatwc_h_ranges_index))
    }
    return(x)
}

# Make all the rasters (using max-stage)
library(parallel)
my_cluster = makeCluster(MC_CORES)
export_job = clusterExport(my_cluster, varlist=ls())
# Combined max-stage and arrival time into one batch of work
combined_files = c(max_stage_rasters_touching_warning_zone, 
    gsub('max_stage', 'arrival_time', max_stage_rasters_touching_warning_zone))
# Spread the jatwc_h_ranges_index over parallel jobs too
combined_indices = seq(1, length(asfm$JATWC_H_ranges))
combined_data = expand.grid(combined_files, combined_indices)
names(combined_data) = c('raster_tile', 'jatwc_h_ranges_index')
# Pack into a list
combined_data_as_list = vector(mode='list', length=nrow(combined_data))
for(i in 1:nrow(combined_data)) combined_data_as_list[[i]] = combined_data[i,]

all_results = parLapplyLB(my_cluster, combined_data_as_list, parallel_fun_combined, chunk.size=1)
stopCluster(my_cluster)
