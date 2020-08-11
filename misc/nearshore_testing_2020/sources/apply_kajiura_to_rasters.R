library(rptha)

recreate_existing_files = TRUE # If TRUE, remake all files. If FALSE, only remake ones that don't exist.

## Elevation DEM (effects kajiura smoothing kernel scale)
elevation_rast = '../elevation/derived_for_model/global/ptha18/merged_gebco_ga250_dem_patched.tif'
#elevation_rast = '/media/gareth/Windows7_OS/Users/gareth/Documents/work/AustPTHA/MODELS/AustPTHA_c/DATA/ELEV/merged_dem/merged_gebco_ga250_dem_patched.tif'

## Raster files containing Okada deformation, which should have kajiura smoothing applied.
raster_files = c(
    'Chile1960/FujiSatake2013/Fuji_chile1960_sources_SUM.tif',
    # The Chile 1960 Ho-et-al source does not require Kajiura (it was anyway
    # derived with water-surface unit-sources)
    ##### 'Chile1960/HoEtAl2019/Ho_Chile1960_initial_displacement.tif',
    'Sumatra2004/FujiSatake2007/Fuji_andaman2004_unit_sources_SUM.tif',
    'Sumatra2004/LoritoEtAl2010/Lorito_sumatra2004_sources_SUM.tif',
    'Sumatra2004/PiatanesiLorito2007/Piatanesi_sumatra2004_sources_SUM.tif',
    'Chile2010/FujiSatake2013/Fuji_chile2010_sources_SUM.tif',
    'Chile2010/LoritoEtAl2011/Lorito_chile2010_sources_SUM.tif',
    'Tohoku2011/SatakeEtAl2013/Satake_Tohoku11_sources_SUM.tif',
    'Tohoku2011/YamakaziEtAl2018/yamakazi18_Tohoku11_sources_SUM.tif',
    'Tohoku2011/RomanoEtAl2015/Tohoku2011_Romano_source.tif',
    'Chile2015/WilliamsonEtAl2017/Williamson_chile_sources_SUM.tif',
    'Chile2015/RomanoEtAl2016/Illapel_2015_Romano.tif')

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

        writeRaster(smoothed, smoothed_out_file, options=c('COMPRESS=DEFLATE'))
    }
}

