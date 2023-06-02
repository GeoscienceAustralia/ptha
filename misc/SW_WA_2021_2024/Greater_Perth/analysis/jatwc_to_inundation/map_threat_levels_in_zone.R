#
# Make rasters that correspond to JATWC zones. Run after the other script, like:
#
#     Rscript compute_scenario_statistics_in_zone.R 'Perth Coast'
#     Rscript map_threat_levels_in_zone.R 'Perth Coast'
#

library(sf)
library(raster)

#
# INPUTS
#

# ATWS Zones
atws_zones_file = normalizePath('ATWS_ZONES/ATWS_Zones_V2_2_4/ATWS_Zones_V2_2_4.shp')
atws_zones = read_sf(dsn=atws_zones_file, layer=gsub('.shp', '', basename(atws_zones_file), fixed=TRUE))
# Was missing a projection, fix that
st_crs(atws_zones) = st_crs('EPSG:4326')

#ATWS_ZONE_NAME = 'Perth Coast'
ATWS_ZONE_NAME = commandArgs(trailingOnly=TRUE)[1]
warning_zone = atws_zones[which(atws_zones$ATWS_Zones == ATWS_ZONE_NAME),]
if(nrow(warning_zone) != 1) stop('ATWS_ZONE_NAME should match exactly one value in attribute ATWS_Zones')

# Greenslade et al (2020) discuss relation between warning categories and the model wave height 
# statistic (termed 'JATWC_H' in this script):
#     Greenslade, D. J. M.; Uslu, B.; Allen, S. C. R.; Kain, C. L.; Wilson, K. M. &
#     Power, H. E. Evaluation of Australian tsunami warning thresholds using
#     inundation modelling Pure Appl. Geophys., Springer Science and Business Media
#     LLC, 2020, 177, 1425-1436
# 
JATWC_H_ranges = list(land_warning = c(0.55, 999999),
                      marine_warning = c(0.2, 0.55),
                      no_threat = c(-1, 0.2), 
                      major_land_warning = c(1.5, 999999),
                      minor_land_warning = c(0.55, 1.5))

# Find model rasters that touch the warning zone
all_scenario_raster_tars = normalizePath(Sys.glob(
    '../../swals/OUTPUTS/ptha18-GreaterPerth-sealevel60cm-reviseddomain-highres/*/*/raster_output_files.tar'))

# How many cores to use?
MC_CORES = 48

#
# END INPUTS
#

# Determine which rasters from the multidomain are touching the warning zone
get_rasters_touching_warning_zone<-function(
    scenario_raster_tar, 
    warning_zone, 
    STARTING_DIR,
    skip_if_cellsize_above_threshold=0.9/(60*9), # For the greater perth model, this excludes the coarsest 2 domain levels
    warning_zone_buffer_degrees=0.1){

    # Get filenames without extraction for speed
    all_raster_basenames = system(paste0('tar --list -f ', scenario_raster_tar), intern=TRUE)
    all_tifs = paste0('/vsitar/', scenario_raster_tar, '/', all_raster_basenames[grep('max_stage_domain_', all_raster_basenames)])

    overlapping = rep(NA, length(all_tifs))

    # Generous buffer in warning zone to include tifs
    warning_zone_buffer = st_buffer(warning_zone, dist=warning_zone_buffer_degrees) # Buffer dist in degrees

    # For each raster, figure out if the extent touches the warning zone
    #for(i in 1:length(all_tifs)){
    parfun<-function(i){
        r1_tiff = raster(all_tifs[i])
        r1_dx = res(r1_tiff)

        if(any(r1_dx > skip_if_cellsize_above_threshold) ){
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
    JATWC_H_range = c(0.5, 999999),
    operation='maxima'){

    # For each scenarios, domain rasters are stored in tar archives -- GDAL can read inside tar archives using the right
    # file name notation ('/vsitar/')
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
        # FIXME: Consider stars -- I think like `t(read_stars(all_raster_files[i])[[1]])`
        r1 = as.matrix(raster(all_raster_files[i]))
        stopifnot(all(dim(r1) == dim(raster_extrema_mat)))

        # Fix for imperfectly recognized missing data values
        r1[r1 < -3e+38] = NA
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

    # FIXME: Consider stars::st_as_stars
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
rasters_touching_warning_zone = get_rasters_touching_warning_zone(
    all_scenario_raster_tars[1], warning_zone, STARTING_DIR)

# Make rasters for all JATWC_H_ranges on a single tile
tile_function<-function(raster_tile, operation){

    for(i in 1:length(JATWC_H_ranges)){
        domain_warning = get_raster_extremes_stratified_by_H(
            raster_tile, all_scenario_JATWC_H, JATWC_H_ranges[[i]], operation=operation)
        output_raster = paste0(names(JATWC_H_ranges)[i], '_', raster_tile)
        writeRaster(domain_warning, output_raster, options=c('COMPRESS=DEFLATE'), overwrite=TRUE)
    }

    rm(domain_warning, output_raster); gc()
    return(0)
}

# Wrap in 'try' to stop failure killing all parallel calculations
parallel_fun_max<-function(raster_tile){
    x = try(tile_function(raster_tile, operation='maxima'))
    return(x)
}

# Make all the rasters (using max-stage)
library(parallel)
all_results = mclapply(rasters_touching_warning_zone, parallel_fun_max, mc.cores=MC_CORES)

# Find the min arrival time
parallel_fun_min<-function(raster_tile){
    x = try(tile_function(raster_tile, operation='minima'))
    return(x)
}
all_results = mclapply(gsub('max_stage', 'arrival_time', rasters_touching_warning_zone), parallel_fun_min, mc.cores=MC_CORES)


