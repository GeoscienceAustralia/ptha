#
# Quick calculation of flood hazard categories
#
source('flood_hazard_categories_FB03.R')
input_par = commandArgs(trailingOnly=TRUE)[1]

if(input_par == "") input_par = '1in2500_84pc'

if(input_par == "1in2500_84pc"){
    # Deliberately use depth-above-initial-condition to "kind-of" reduce issues with too-high initial stage.
    string_matching_depth_files = "nsw_full_coast_MIS_max_depth_above_initial_condition_1in2500_84pc_at_sites_with_elevation_above_0/*.tif" 
    string_matching_speed_files = "nsw_full_coast_MIS_max_speed_1in2500_84pc/*.tif" 
    string_matching_flux_files = "nsw_full_coast_MIS_max_flux_1in2500_84pc/*.tif" 
    output_dir = 'nsw_full_coast_MIS_flood_hazard_categories_1in2500_84pc'

}else if(input_par == "1in250_50pc"){

    string_matching_depth_files = "nsw_full_coast_MIS_max_depth_above_initial_condition_1in250_50pc_at_sites_with_elevation_above_0/*.tif" 
    string_matching_speed_files = "nsw_full_coast_MIS_max_speed_1in250_50pc/*.tif" 
    string_matching_flux_files = "nsw_full_coast_MIS_max_flux_1in250_50pc/*.tif" 
    output_dir = 'nsw_full_coast_MIS_flood_hazard_categories_1in250_50pc'

}else{
    stop(paste0('unknown input_par: ', input_par))
}

depth_files = Sys.glob(string_matching_depth_files)
vd_files = Sys.glob(string_matching_flux_files)
speed_files = Sys.glob(string_matching_speed_files)
#output_dir = 'trial_hazard_categories'

ldf = length(depth_files)
stopifnot((ldf > 0) & (length(vd_files) == ldf) & (length(speed_files) == ldf) )

dir.create(output_dir, showWarnings=FALSE)

# Convert the data into a convenient structure for parallel execution
parallel_data = mapply(function(vd, depth, speed) list(vd=vd, depth=depth, speed=speed), 
    vd=vd_files, depth=depth_files, speed=speed_files, SIMPLIFY=FALSE, USE.NAMES=FALSE)

library(parallel)
MC_CORES=16
mclapply(parallel_data, function(x){
    vd = rast(x$vd)
    depth = rast(x$depth)
    speed = rast(x$speed)
    hazard_categories = compute_flood_hazard_categories(vd, depth, speed)

    # Skip if we have missing data for any input variable.
    hazard_categories[is.na(depth) | is.na(vd) | is.na(speed)] = NA

    output_file = paste0(output_dir, '/', gsub("_max_speed_", "_hazard_categories_", basename(x$speed)))
    writeRaster(hazard_categories, file=output_file, gdal=c('COMPRESS=DEFLATE'), overwrite=TRUE)
    return(0)
    }, 
    mc.cores=MC_CORES,  mc.preschedule=TRUE)
