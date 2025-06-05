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
asfm$SCENARIO_AMBIENT_SEA_LEVEL = 0.0 # m
asfm$WETTOL = 0.01 # m
asfm$jatwc_depth_limit = 20 # m
asfm$jatwc_dlon_dlat = 4/60 # degrees
asfm$scenario_elevation_raster_file = '../OUTPUTS/tohoku2011-full-ambient_sea_level_0/RUN_20241120_163301922/elevation0.vrt'
asfm$atws_zones = st_read('../../analysis/jatwc_to_inundation/ATWS_ZONES/ATWS_Zones_V2_2_4/ATWS_Zones_V2_2_4.shp')
asfm$all_scenario_raster_dirs = c(
    '../OUTPUTS/solomon2007_1_19-full-ambient_sea_level_0/RUN_20241120_163301863',
    '../OUTPUTS/tohoku2011-full-ambient_sea_level_0/RUN_20241120_163301922',
    '../OUTPUTS/chile2010_fuji-full-ambient_sea_level_0/RUN_20250108_134230646'
)

# check that all the directories exist and have max_stage*.tif files
for (scenario_raster_dir in asfm$all_scenario_raster_dirs) {
    if (!dir.exists(scenario_raster_dir)) {
        stop(paste('Directory does not exist:', scenario_raster_dir))
    }
    if (length(Sys.glob(file.path(scenario_raster_dir, 'max_stage*.tif'))) == 0) {
        stop(paste('No max_stage*.tif files in directory:', scenario_raster_dir))
    }
}


# The ATWS warning zone we work on (must match value in ATWS_Zones column of ATWS shapefile)
ATWS_ZONE_NAME = 'Capricornia Coast'

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

#
# END INPUTS
#

#
# The JATWC warnings are based on '95th percentile stage within the warning zone'. We compute this by selecting points
# in the warning zone (roughly matching JATWC model resolution and depth constraints) and then extracting the max-stage.
# This function creates the points, using caching so that expensive computations only happen once.
#
get_offshore_points_to_compute_JATWC_H<-function(warning_zone, scenario_elevation_file, 
    depth_limit, dlon_dlat, jatwc_dir='.'){

    # Cache points to this file
    output_file = paste0(
        jatwc_dir,
        'comparison-points-for-JATWC-H_', 
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
# Creation of vrt, calculation of JATWC_H
#
make_max_stage_rast_and_get_JATWC_H<-function(scenario_raster_dir, sample_points){
    # Make a VRT with the max-stage data
    tmp_vrt = './tmp_all_max_stage_.vrt'
    cmd = paste0('gdalbuildvrt -resolution highest ', tmp_vrt, ' ', scenario_raster_dir, '/max_stage*.tif')
    print(cmd)
    system(cmd, ignore.stdout=TRUE)

    # Get the JATWC_H 
    scenario_max_stage = raster(tmp_vrt)
    jatwc_H = get_max_stage_percentile_in_zone(scenario_max_stage, sample_points)

    # event_name from filepath
    if (grepl('solomon2007', scenario_raster_dir)) {
        event_name = 'Solomon2007'
    } else if (grepl('tohoku2011', scenario_raster_dir)) {
        event_name = 'Tohoku2011'
    } else if (grepl ('chile2010', scenario_raster_dir)) {
        event_name = 'chile2010'
    } else {
        stop('Unknown event')
    }
    output = data.frame(model_dir = scenario_raster_dir, jatwc_H = as.numeric(jatwc_H), event_name = event_name)

    # Cleanup
    rm(scenario_max_stage, jatwc_H); gc()
    file.remove(tmp_vrt)

    return(output)
}

#
# Main program here
#

ATWS_Zone_name_nospace = gsub(' ', '-', ATWS_ZONE_NAME)

# Get Points at which we can compute JATWC_H
jatwc_dir = '../../analysis/jatwc_to_inundation/Inundation_zones/Capricornia-Coast'
sample_points = get_offshore_points_to_compute_JATWC_H(
    warning_zone, scenario_elevation_raster_file, asfm$jatwc_depth_limit, asfm$jatwc_dlon_dlat, jatwc_dir=jatwc_dir)
stopifnot(max(sample_points[,3]) <= (-1*asfm$jatwc_depth_limit))

# Compute JATWC_H for all scenarios in a dataframe
all_JATWC_H = data.frame()
for (i in 1:length(asfm$all_scenario_raster_dirs)) {
    scenario_raster_dir = asfm$all_scenario_raster_dirs[i]
    result = make_max_stage_rast_and_get_JATWC_H(scenario_raster_dir, sample_points)
    all_JATWC_H = rbind(all_JATWC_H, result)
}

output_file = paste0('all-JATWC-H_', ATWS_Zone_name_nospace, '.RDS')
saveRDS(all_JATWC_H, file=output_file)
write.csv(all_JATWC_H, paste0('all-JATWC-H_', ATWS_Zone_name_nospace, '.csv'), row.names=FALSE)
