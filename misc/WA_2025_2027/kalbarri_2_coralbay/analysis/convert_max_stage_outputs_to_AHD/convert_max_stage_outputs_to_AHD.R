#
# For the max-stage output products, convert the vertical units from "m above local HAT" to "m above AHD".
#
# The products include:
# - max stage 1/2500 84%
# - marine_warning_max_stage
#

library(terra)

# Tidal adjustment used by SWALS
tidal_adjustment_rasters = Sys.glob('../jatwc_to_inundation/tidal_adjustment/*.tif')

# max stage, 1/2500 84%
max_stage_1in2500_84pc_rasters = Sys.glob('../probabilistic_inundation/kalbarri2coralbay_highres_domains_max_stage_percentile_0.84_exrate_0.0004_hazard/*.tif')
# .... and dir for outputs
max_stage_1in2500_84pc_output_dir = 'kalbarri2coralbay_highres_domains_max_stage_AHD_percentile_0.84_exrate_0.0004_hazard'
nothing = dir.create(max_stage_1in2500_84pc_output_dir, showWarnings=FALSE)

# 'marine warning max stage' for each ATWS zone
available_inundation_zones = c('Gascoyne-Coast', 'Geraldton-Coast', 'Ningaloo-Coast')
marine_warning_max_stage_rasters = lapply(available_inundation_zones,
    function(x) Sys.glob(paste0('../jatwc_to_inundation/Inundation_zones/', x, '/marine_warning_max_stage_domain_*.tif')))
names(marine_warning_max_stage_rasters) = available_inundation_zones
# ... and dir for outputs
marine_warning_max_stage_output_dirs = lapply(available_inundation_zones, 
    function(x){paste0('marine_warning_max_stage_AHD/', x, '/')})
names(marine_warning_max_stage_output_dirs) = available_inundation_zones
nothing = lapply(marine_warning_max_stage_output_dirs, function(x) dir.create(x, showWarnings=FALSE, recursive=TRUE))

#
# The code loops over domains, applying the tidal adjustment on the corresponding
# domain. So we need a way to find the domain index.
#
find_domain_index_from_filename=function(my_raster_file){

    stopifnot(length(my_raster_file) == 1) # Not vectorised

    # Convenience function for getting the index
    extract_ind=function(splitter) as.numeric(gsub('.tif', '', strsplit(my_raster_file, splitter)[[1]][2]), fixed=TRUE) 

    if(grepl('_domain_index_', my_raster_file)){
        # Some raster filenames end in _domain_index_1234.tif where 1234 is the index we want
        domain_index = extract_ind('_domain_index_')
    }else if(grepl('_domain_', my_raster_file)){
        # Some raster filenames end in _domain_1234.tif where 1234 is the index we want
        domain_index = extract_ind('_domain_')
    }else{
        stop('Could not identify domain index')
    }
    domain_index
}

# Store tidal adjustment indices which we will later match within a loop
tidal_adjustment_domain_inds = sapply(tidal_adjustment_rasters, find_domain_index_from_filename)

#
# Adjust the max-stage 1in2500 rasters
#
adjust_max_stage_1in2500_raster_i<-function(i){

    r1_file = max_stage_1in2500_84pc_rasters[i]

    # Find the corresponding tidal adjustment tif
    r1_ind = find_domain_index_from_filename(r1_file)
    tidal_ind = which(r1_ind == tidal_adjustment_domain_inds)
    if(length(tidal_ind) != 1) stop(paste0('No unique match for ', r1_file))

    r1 = rast(r1_file) # Result in 'm above max tide'
    t1 = rast(tidal_adjustment_rasters[tidal_ind]) # 'max tide above AHD'

    # Convert max stage to units of AHD
    r1_AHD = r1 + t1

    # Save it
    r1_AHD_filename = paste0(max_stage_1in2500_84pc_output_dir, '/', gsub('_max_stage_', '_max_stage_AHD_', basename(r1_file)))
    writeRaster(r1_AHD, file=r1_AHD_filename, gdal=c('COMPRESS=DEFLATE'), overwrite=TRUE)
}
for(i in 1:length(max_stage_1in2500_84pc_rasters)) adjust_max_stage_1in2500_raster_i(i)

#
# Adjust the marine warning max stage rasters
#
adjust_marine_warning_max_stage_raster_i<-function(i, inundation_zone){

    r1_file = marine_warning_max_stage_rasters[[inundation_zone]][i]

    # Find the corresponding tidal adjustment tif
    r1_ind = find_domain_index_from_filename(r1_file)
    tidal_ind = which(r1_ind == tidal_adjustment_domain_inds)
    if(length(tidal_ind) != 1) stop(paste0('No unique match for ', r1_file))

    r1 = rast(r1_file) # Result in 'm above max tide'
    t1 = rast(tidal_adjustment_rasters[tidal_ind]) # 'max tide above AHD'

    # Convert max stage to units of AHD
    r1_AHD = r1 + t1

    # Save it
    r1_AHD_filename = paste0(marine_warning_max_stage_output_dirs[[inundation_zone]], '/', gsub('max_stage_', 'max_stage_AHD_', basename(r1_file)))
    writeRaster(r1_AHD, file=r1_AHD_filename, gdal=c('COMPRESS=DEFLATE'), overwrite=TRUE)
}
for(inundation_zone in available_inundation_zones){
    N = length(marine_warning_max_stage_rasters[[inundation_zone]])
    for(i in 1:N) adjust_marine_warning_max_stage_raster_i(i, inundation_zone) 
}
