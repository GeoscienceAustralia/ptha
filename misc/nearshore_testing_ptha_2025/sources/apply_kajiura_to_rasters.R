library(rptha)

recreate_existing_files = TRUE # If TRUE, remake all files. If FALSE, only remake ones that don't exist.

## Elevation DEM (effects kajiura smoothing kernel scale)
elevation_rast = '../elevation/derived_for_model/global/ptha18/merged_gebco_ga250_dem_patched.tif'
#elevation_rast = '/media/gareth/Windows7_OS/Users/gareth/Documents/work/AustPTHA/MODELS/AustPTHA_c/DATA/ELEV/merged_dem/merged_gebco_ga250_dem_patched.tif'

## Raster files containing Okada deformation, which should have kajiura smoothing applied.
raster_files = c(
    # The Chile 1960 Ho-et-al source does not require Kajiura (it was anyway
    # derived with water-surface unit-sources)
    ##### 'Chile1960/HoEtAl2019/Ho_Chile1960_initial_displacement.tif',
    # Sumatra 2004 is in the time-varying sources (separate Kajiura code therein)
    # Sumatra 2005 was dealt with separately (inherited from WA tsunami project, but uses a per-unit-source Kajiura code like the time-varying sources)
    # Java 2006 was in the time-varying sources (separate Kajiura code therein)
    'Solomon2007/Wei_2015/Wei_S2_Solomon2007_source_SUM.tif',
    'Sumatra2007/FujiiSatake2007/FujiiSatake08_2007_Sumatra_SUM.tif',
    'Puysegur2009/Bevan2010/Bevan_Puysegur2009_source_SUM.tif',
    'Chile2010/LoritoEtAl2011/Lorito_chile2010_sources_SUM.tif',
    'Tohoku2011/YamazakiEtAl2018Fixed/yamazaki18_Tohoku11_sources_SUM.tif',
    'Chile2014/An_2014/An_Chile2014_source_SUM.tif',
    'Chile2015/WilliamsonEtAl2017/Williamson_chile_sources_SUM.tif',
    'NewHebrides2021/GusmanEtAl/Gusman_NewHebrides_sources_SUM.tif',
    # The Kermadec2021 Romano source model already has Kajiura in the unit sources, do not apply again.
    #####'KermadecTonga2021/Romano2021/Kermadec2021_Romano.tif
    'Sandwich2021/Roger2024_GCMT_source/RogerGCMT_2024_Sandwich2021_SUM.tif'
)

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

