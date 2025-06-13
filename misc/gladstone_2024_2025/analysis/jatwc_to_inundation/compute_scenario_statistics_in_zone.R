#
# Compute the JATWC parameter for all scenarios.
#
# The JATWC parameter 'H' is the 95th percentile of tsunami maxim in a zone on
# the T2 model grid, where model depths are > 20m, and the grid size is 4
# arcmin.
#
# Run like:
#     Rscript compute_scenario_statistics_in_zone.R 'Perth Coast'
#
#
library(sf)
library(raster)

#
# INPUTS
#
asfm = new.env()
source('application_specific_metadata.R', local=asfm)

# The ATWS warning zone we work on (must match value in ATWS_Zones column of ATWS shapefile)
ATWS_ZONE_NAME = commandArgs(trailingOnly=TRUE)[1] # 'Perth Coast'

##
## END INPUTS
##

# A reference scenario used to get elevation for all models
#scenario_elevation = raster(asfm$scenario_elevation_raster_file)
scenario_elevation_raster_file = normalizePath(asfm$scenario_elevation_raster_file)

# ATWS Zones Shapefile
atws_zones = asfm$atws_zones
warning_zone = atws_zones[which(atws_zones$ATWS_Zones == ATWS_ZONE_NAME),]
if(nrow(warning_zone) != 1) stop('ATWS_ZONE_NAME should match exactly one value in attribute ATWS_Zones')

# How many cores for parallel computation?
MC_CORES = asfm$DEFAULT_MC_CORES

#
# END INPUTS
#

#
# The JATWC warnings are based on '95th percentile stage within the warning zone'. We compute this by selecting points
# in the warning zone (roughly matching JATWC model resolution and depth constraints) and then extracting the max-stage.
# This function creates the points, using caching so that expensive computations only happen once.
#
get_offshore_points_to_compute_JATWC_H<-function(warning_zone, scenario_elevation_file, 
    depth_limit, dlon_dlat){

    # Cache points to this file
    output_file = paste0('comparison-points-for-JATWC-H_', 
        warning_zone$ATWS_Zones, '_', 
        depth_limit, '_', 
        signif(dlon_dlat, 5), '.RDS')
    output_file = gsub(' ', '-', output_file)

    if(file.exists(output_file)){
        # Quick version
        sample_points = readRDS(output_file)
    }else{
        # Find points in the warning zone that are not too deep

        # Firstly make a grid of points
        wz_bbox = st_bbox(warning_zone)
        slon = seq(floor(wz_bbox$xmin), ceiling(wz_bbox$xmax), by=dlon_dlat)
        slat = seq(floor(wz_bbox$ymin), ceiling(wz_bbox$ymax), by=dlon_dlat)
        # Candidate grid points
        grid_pts = expand.grid(slon, slat)

        # Get their elevation -- consider stars::st_extract
        scenario_elevation = raster(scenario_elevation_file)
        grid_pts_elev = extract(scenario_elevation, grid_pts)

        dim(grid_pts_elev) = c(length(slon), length(slat))
        grid_elev = raster(t(grid_pts_elev), 
            xmn =min(slon), xmx=max(slon), 
            ymn=min(slat), ymx=max(slat), crs=proj4string(scenario_elevation))
        grid_elev = flip(grid_elev, direction='y')
        # Consider stars::st_as_stars
        candidate_points = rasterToPoints(grid_elev)

        # For the general case we might want to remove isolated clumps of points here
        # (which would not be connected in the JATWC model). 
        # We can check for this in the output PNG file -- to date it hasn't been an issue.

        # Find raster cells in polygon
        wz_rast = rasterize(warning_zone, grid_elev)
        is_in_poly = extract(wz_rast, candidate_points[,1:2])

        # Keep points in the polygon that are deep enough. Elevation (m) is in
        # candidate_points[,3], ocean values are negative.
        keep = which((candidate_points[,3] < (-1*depth_limit)) &
                     (!is.na(is_in_poly)))
        
        # Get lon,lat,elev at points we want to keep
        sample_points = candidate_points[keep,]

        # Save to a file for reading later
        saveRDS(sample_points, output_file)

        # Make a PNG so we can see where the points are
        library(sp) 
        png('Points_for_JATWC_H.png', width=6, height=6, units='in', res=300)
        plot(as_Spatial(warning_zone), col='green', 
             main = paste0(warning_zone$ATWS_Zones, 
                 ': Points used to estimate JATWC H \n (95th percentile tsunami maxima with depth > 20m)'), 
             asp=1, axes=TRUE, cex.main=0.9)
        plot(as_Spatial(atws_zones), border='black', col=NA, add=TRUE)
        points(sample_points[,1:2], pch=19, cex=0.2, col='red')
        dev.off()
    }

    to_remove = setdiff(ls(), 'sample_points')
    rm(list=to_remove); gc()

    return(sample_points)
}

#
# Get the 95th percentile max-stage at sample points in the warning zone,
# excluding dry sites and locations that the tsunami never reaches.
#
get_max_stage_percentile_in_zone<-function(scenario_max_stage, sample_points, prob=0.95){

    # Extract stage/elevation at the sampled points
    sample_stage = extract(scenario_max_stage, sample_points[,1:2])
    sample_elev = sample_points[,3]

    # Exclude points with NA stage, or dry points, or points that never exceed their initial stage
    keep = which( (sample_stage > (sample_elev + asfm$WETTOL)) &
                  (sample_stage > (asfm$SCENARIO_AMBIENT_SEA_LEVEL + asfm$WETTOL)) &
                  !is.na(sample_stage))

    # JATWC H parameter
    jatwc_H_parameter = quantile(sample_stage[keep], probs=prob, type=6) - asfm$SCENARIO_AMBIENT_SEA_LEVEL

    rm(sample_stage, sample_elev, keep); gc()

    return(c(jatwc_H_parameter))
}

#
# Wrapper that deals with untarring of raster files, creation of vrt, calculation of JATWC_H, and cleanup
#
make_max_stage_rast_and_get_JATWC_H<-function(scenario_raster_tar, sample_points, STARTING_DIR, extra_dir_tag){

    # Ensure we always arrive back at STARTING_DIR
    on.exit(setwd(STARTING_DIR))

    # Use a local temporary directory. Append extra_dir_tag in case we have
    # multiple tars with the same parent directory name (which can happen with e.g. 
    # multiple batches of scenarios using links for repeated scenarios)
    local_dir = paste0(basename(dirname(scenario_raster_tar)), '_', extra_dir_tag)

    dir.create(local_dir)
    setwd(local_dir)

    ##
    # FIXME: For the code below, in future it may be faster to avoid creating
    # a VRT, and just do the extraction at specific rasters, read from inside
    # their tar archive.
    ##

    # Extract the max-stage domain tifs to local_dir
    system(paste0('tar -xf ', scenario_raster_tar, ' --wildcards "max_stage_domain*.tif"'))

    # Make a VRT with the max-stage data
    # Consider stars:: -- read_stars(st_mosaic( list-of-files, options=c('-resolution', 'highest')))
    system('gdalbuildvrt -resolution highest all_max_stage.vrt max_stage*.tif', ignore.stdout=TRUE)

    # Get the JATWC_H 
    scenario_max_stage = raster('all_max_stage.vrt')
    jatwc_H = get_max_stage_percentile_in_zone(scenario_max_stage, sample_points)

    output = data.frame(model_dir = dirname(scenario_raster_tar), jatwc_H = as.numeric(jatwc_H))

    setwd(STARTING_DIR)

    # Remove the files
    unlink(local_dir, recursive=TRUE)

    rm(scenario_max_stage, jatwc_H); gc()

    return(output)
}

#
# Main program here
#

ATWS_Zone_name_nospace = gsub(' ', '-', ATWS_ZONE_NAME)

# Make space for outputs
working_dir = paste0('Inundation_zones/', ATWS_Zone_name_nospace)
dir.create(working_dir, showWarnings=FALSE, recursive=TRUE)
setwd(working_dir)

# Useful to return here
STARTING_DIR = getwd()

# Get Points at which we can compute JATWC_H
sample_points = get_offshore_points_to_compute_JATWC_H(
    warning_zone, scenario_elevation_raster_file, asfm$jatwc_depth_limit, asfm$jatwc_dlon_dlat)
stopifnot(max(sample_points[,3]) <= (-1*asfm$jatwc_depth_limit))

# Wrapper to run in parallel
parallel_fun<-function(i){
    library(sf)
    library(raster)
    scenario_raster_tar = asfm$all_scenario_raster_tars[i]
    result = try(make_max_stage_rast_and_get_JATWC_H(scenario_raster_tar, sample_points, STARTING_DIR, extra_dir_tag = i))
    return(result) 
}

library(parallel)
my_cluster = makeCluster(MC_CORES)
export_data = clusterExport(my_cluster, ls())
#all_JATWC_H = mclapply(1:length(asfm$all_scenario_raster_tars), parallel_fun, mc.cores=MC_CORES)
all_JATWC_H = parLapplyLB(my_cluster, 1:length(asfm$all_scenario_raster_tars), parallel_fun, chunk.size=1)
stopCluster(my_cluster)
output_file = paste0('all-JATWC-H_', ATWS_Zone_name_nospace, '.RDS')
saveRDS(all_JATWC_H, file=output_file)
