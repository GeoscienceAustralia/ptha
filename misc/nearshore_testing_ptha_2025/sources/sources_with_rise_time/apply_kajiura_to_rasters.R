# Apply Kajiura smoothing to some rasters
#
# Can run in parallel by setting off many jobs like
#
# Rscript apply_kajiura_to_rasters 1 3
# Rscript apply_kajiura_to_rasters 2 3
# Rscript apply_kajiura_to_rasters 3 3

library(rptha)

recreate_existing_files = FALSE # If TRUE, remake all files. If FALSE, only remake ones that don't exist.

## Elevation DEM (effects kajiura smoothing kernel scale)
#elevation_rast = '../../../../../MODELS/AustPTHA_c/DATA/ELEV/merged_dem/merged_gebco_ga250_dem_patched.tif'
elevation_rast = '/media/gareth/Windows7_OS/Users/gareth/Documents/work/AustPTHA/MODELS/AustPTHA_c/DATA/ELEV/merged_dem/merged_gebco_ga250_dem_patched.tif'

raster_files = c(
    #Sys.glob('Sumatra2004/FujiSatake2007/Fuji_andaman2004_*.tif'),
    Sys.glob('Java2006/FujiiSatake06/FujiiSatake06_*.tif')
)
# Avoid applying smoothing to files that we already smoothed
k = grep('KAJIURA_SMOOTHED', raster_files, fixed=TRUE)
if(length(k) > 0) raster_files = raster_files[-k]

library(parallel)
input_args = as.numeric(commandArgs(trailingOnly=TRUE))
my_inds = splitIndices(length(raster_files), input_args[2])[[input_args[1]]]
raster_files = raster_files[my_inds]


all_file_paths = c(elevation_rast, raster_files)
file_paths_exist = file.exists(all_file_paths)
if(!all(file_paths_exist)){
    print(paste0(file_paths_exist, ' ', all_file_paths))
    print('Could not find some files -- see files corresponding to "FALSE" above')
}


for(i in 1:length(raster_files)){
    initial = raster(raster_files[i])
    
    smoothed_out_file = paste0(gsub('.tif', '', raster_files[i], fixed=TRUE), '_KAJIURA_SMOOTHED.tif')
    if(recreate_existing_files | !file.exists(smoothed_out_file)){

        # For smoothing we need to linearize the coordinates about some origin
        # Choose the middle of the input raster. 
        extent_initial = extent(initial)
        new_origin = 0.5*c((extent_initial@xmin + extent_initial@xmax), 
                           (extent_initial@ymin + extent_initial@ymax))

        # Do the smoothing
        smoothed = kajiura_smooth_raster(initial, new_origin = new_origin, 
            elevation_raster = elevation_rast, spherical_input=TRUE)

        writeRaster(smoothed, smoothed_out_file, options=c('COMPRESS=DEFLATE'), overwrite=TRUE)
    }
}
