# Make the SWALS elevation file preference list.
# The original list was created on my home machine (interactively ordering files in GIS). 
# This script reads that list, then updates the file paths to point to the locations on NCI

files = read.csv('ordered_files_home_machine_rev1/Elevation_rasters_scraped_from_QGIS_file.csv')

new_files = files$files

# New locations for the "Gap_filling_estuaries" folder.
new_files = gsub(
    '/media/gareth/Windows7_OS/Users/gareth/Documents/work/Inundation_tsunami/bathymetry_merging/NSW/Gap_filling_estuaries/',
    '/g/data/w85/tsunami/MODELS/inundation/DATA/ELEV/NSW/Gap_filling_estuaries_2023_08/',
    new_files, fixed=TRUE)

# New location for first priority manual patches
new_files = gsub(
    '/media/gareth/Windows7_OS/Users/gareth/Documents/work/Inundation_tsunami/bathymetry_merging/NSW/patches_with_first_priority/',
    '/g/data/w85/tsunami/MODELS/inundation/DATA/ELEV/NSW/extra_manual_patches/patches_with_first_priority/',
    new_files, fixed=TRUE)

# New location for the coastal lidar with patching
new_files = gsub(
    '/media/gareth/Windows7_OS/Users/gareth/Documents/work/Inundation_tsunami/bathymetry_merging/NSW/coastal_lidar_masked_by_WoFS_patched/',
    '/g/data/w85/tsunami/MODELS/inundation/DATA/ELEV/NSW/coastal_lidar_masked_by_WoFS_patched/',
    new_files, fixed=TRUE)

# New location for patched 2018 topobathy
new_files = gsub(
    '/media/gareth/Windows7_OS/Users/gareth/Documents/work/Inundation_tsunami/bathymetry_merging/NSW/bathy2018_patched/',
    '/g/data/w85/tsunami/MODELS/inundation/DATA/ELEV/NSW/bathy2018_patched/',
    new_files, fixed=TRUE)

# New location for the patched 2018 topobathy with holes filled.
new_files = gsub(
    '/media/gareth/Windows7_OS/Users/gareth/Documents/work/Inundation_tsunami/bathymetry_merging/NSW/bathy2018_patched_fillnodata/',
    '/g/data/w85/tsunami/MODELS/inundation/DATA/ELEV/NSW/bathy2018_patched_fillnodata/',
    new_files, fixed=TRUE)

# New location for batemans merged highres
new_files = gsub(
    '/media/gareth/Windows7_OS/Users/gareth/Documents/work/Inundation_tsunami/bathymetry_merging/NSW/batemans_merge_highres/',
    '/g/data/w85/tsunami/MODELS/inundation/DATA/ELEV/NSW/batemans_merge_highres/',
    new_files, fixed=TRUE)

# New locations for Wilson and power DEMs
## Botany
new_files = gsub(
    '/media/gareth/Windows7_OS/Users/gareth/Documents/work/My Documents/Documents/Gareth/projects/tsunami/references/DATA/Wilson_and_Power_Sydney_Bathymetry/Botany/botany_0.0001_gcs.txt.tif',
    '/g/data/w85/tsunami/DATA/ELEVATION/NSW/SydneyRegion/botany_0.0001_gcs.txt.tif',
    new_files, fixed=TRUE)
## Sydney
new_files = gsub(
    '/media/gareth/Windows7_OS/Users/gareth/Documents/work/My Documents/Documents/Gareth/projects/tsunami/references/DATA/Wilson_and_Power_Sydney_Bathymetry/Sydney/syd_0.0001_gcs.txt.tif',
    '/g/data/w85/tsunami/DATA/ELEVATION/NSW/SydneyRegion/syd_0.0001_gcs.txt.tif',
    new_files, fixed=TRUE)
## Hawkesbury
new_files = gsub(
    '/media/gareth/Windows7_OS/Users/gareth/Documents/work/My Documents/Documents/Gareth/projects/tsunami/references/DATA/Wilson_and_Power_Sydney_Bathymetry/hawkesbury_revised/hawkesbury_0.0005_gcs.txt.tif',
    '/g/data/w85/tsunami/DATA/ELEVATION/NSW/SydneyRegion/hawkesbury_0.0005_gcs.txt.tif',
    new_files, fixed=TRUE)

# New locations for Gold Coast high res DEM.
# This has a different md5sum compared to home machine (written with different library versions?).
new_files = gsub(
    '/media/gareth/Windows7_OS/Users/gareth/Documents/work/AustPTHA/DATA/ELEV/QLD/Gold_Coast_DTM/GCCC_2m_DTM_Sept_2020_WGS84_noboundaryzeros_add_channel.tif',
    '/g/data/w85/tsunami/DATA/ELEVATION/Gold_Coast_DTM/GCCC_2m_DTM_Sept_2020_WGS84_noboundaryzeros_add_channel.tif',
    new_files, fixed=TRUE)

# Hunter River good patch DEM
new_files = gsub(
    '/media/gareth/Windows7_OS/Users/gareth/Documents/work/Inundation_tsunami/bathymetry_merging/NSW/good_patch_dems/hunter_river_gridded_WGS84.tif',
    '/g/data/w85/tsunami/MODELS/inundation/DATA/ELEV/NSW/Gap_filling_estuaries_2023_08/newcastle_harbour/upper_channels/create_patch_dem_newcastle/hunter_river_gridded_WGS84.tif',
    new_files, fixed=TRUE)

# Wonboyn lake good patch DEM
new_files = gsub(
    '/media/gareth/Windows7_OS/Users/gareth/Documents/work/Inundation_tsunami/bathymetry_merging/NSW/Gap_filling_estuaries/wonboyn_lake/wonboyn_wgs84.tif',
    '/g/data/w85/tsunami/MODELS/inundation/DATA/ELEV/NSW/Gap_filling_estuaries_2023_08/wonboyn_lake/wonboyn_wgs84.tif',
    new_files, fixed=TRUE)

# Clarence river good patch DEMs
new_files = gsub(
    '/media/gareth/Windows7_OS/Users/gareth/Documents/work/Inundation_tsunami/bathymetry_merging/NSW/Gap_filling_estuaries/clarence_river',
    '/g/data/w85/tsunami/MODELS/inundation/DATA/ELEV/NSW/Gap_filling_estuaries_2023_08/clarence_river',
    new_files, fixed=TRUE)

# Richmond river good patch DEMs
new_files = gsub(
    '/media/gareth/Windows7_OS/Users/gareth/Documents/work/Inundation_tsunami/bathymetry_merging/NSW/Gap_filling_estuaries/richmond_river',
    '/g/data/w85/tsunami/MODELS/inundation/DATA/ELEV/NSW/Gap_filling_estuaries_2023_08/richmond_river',
    new_files, fixed=TRUE)

# Port Kembla, good patch
new_files = gsub(
    '/media/gareth/Windows7_OS/Users/gareth/Documents/work/Inundation_tsunami/bathymetry_merging/NSW/good_patch_dems/PortKembla_high_res.tif',
    '/g/data/w85/tsunami/MODELS/inundation/DATA/ELEV/NSW/highres_coastal_merge/PortKembla_high_res.tif',
    new_files, fixed=TRUE)

# Lord Howe Island, various DEMs
new_files = gsub(
    '/media/gareth/Data2/DATA/Lord_Howe_Island',
    '/g/data/w85/tsunami/DATA/ELEVATION/Lord_Howe_Island',
    new_files, fixed=TRUE)

# Norfolk Island, various DEMS
new_files = gsub(
    "/media/gareth/Data2/DATA/Norfolk_Island/Nearshore_LADS/Derived_product_which_can_be_distributed/",
    "/g/data/w85/tsunami/DATA/ELEVATION/Norfolk_Island/AHO_DEM_derived_from_nearshore_lads_can_be_distributed/", 
    new_files, fixed=TRUE)
new_files = gsub(
    "/media/gareth/Data2/DATA/Norfolk_Island/Norfolk_Island_Lidar/",
    "/g/data/w85/tsunami/DATA/ELEVATION/Norfolk_Island/Norfolk_Island_Lidar/",
    new_files, fixed=TRUE)
new_files = gsub(
   "/media/gareth/Data2/DATA/Norfolk_Island/Norfolk_Island_Nearshore_Coastal_2021/", 
   "/g/data/w85/tsunami/DATA/ELEVATION/Norfolk_Island/Norfolk_Island_Nearshore_Coastal_2021/",
    new_files, fixed=TRUE)
    

# Myall River, good patch -- this is taken care of in other "Gap_filling_estuaries" updates

# Williams river Newcastle -- this is taken care of in other "Gap filling estuaries" updates

# Big victorian DEM
new_files = gsub(
    '/media/gareth/Windows7_OS/Users/gareth/Documents/work/AustPTHA/DATA/ELEV/Vic/',
    '/g/data/w85/tsunami/DATA/ELEVATION/Victoria/Vic/',
    new_files, fixed=TRUE)

# masked_by_WOFS_patched lidar -- taken care of above

# Many "Gap_filling_estuaries" patches are taken care of above.

# Hawkesbury Wilson and Power - taken care of above.

# Great Barrier Reef 30m, 2020 DEM
new_files = gsub(
    '/media/gareth/Data2/DATA/Great_Barrier_Reef_30m_2020/',
    '/g/data/w85/tsunami/DATA/ELEVATION/GreatBarrierReef_30m/Great_Barrier_Reef_30m_2020/',
    new_files, fixed=TRUE)

# Patched 2018 topography with holes filled, treated above

# Transition DEMs
new_files = gsub(
    '/media/gareth/Windows7_OS/Users/gareth/Documents/work/Inundation_tsunami/bathymetry_merging/NSW/ocean_smooth_transition/',
    '/g/data/w85/tsunami/MODELS/inundation/DATA/ELEV/NSW/ocean_smooth_transition/',
    new_files, fixed=TRUE)
new_files = gsub(
    '/media/gareth/Windows7_OS/Users/gareth/Documents/work/Inundation_tsunami/bathymetry_merging/NSW/ocean_smooth_transition_east_coast/',
    '/g/data/w85/tsunami/MODELS/inundation/DATA/ELEV/NSW/ocean_smooth_transition_east_coast/',
    new_files, fixed=TRUE)
new_files = gsub(
    '/media/gareth/Windows7_OS/Users/gareth/Documents/work/Inundation_tsunami/bathymetry_merging/Victoria/',
    '/g/data/w85/tsunami/MODELS/inundation/DATA/ELEV/NSW/ocean_smooth_transition_victoria/Victoria/',
    new_files, fixed=TRUE)

# Patches of the 2023 GA250 DEM.
# Not using this everywhere because
#   A) PTHA18 used the global DEM below -- and initial NSW model deveopment did too -- so results will be perturbed less by sticking to that.
#   B) The 2023 DEM has some artefacts (fake islands) that I'd rather not clean up at this point
new_files = gsub(
    "/media/gareth/Data2/DATA/Australian_Bathymetry_and_Topography_Grid_2023/clip_subsets/",
    "/g/data/w85/tsunami/DATA/ELEVATION/Australia/Australian_Bathymetry_and_Topography_2023_250m/clip_subsets/",
    new_files, fixed=TRUE)

# PTHA18 Global DEM
new_files = gsub(
    '/media/gareth/Windows7_OS/Users/gareth/Documents/work/AustPTHA/MODELS/AustPTHA_c/DATA/ELEV/merged_dem/merged_gebco_ga250_dem_patched.tif',
    '/g/data/w85/tsunami/MODELS/inundation/australia_wide/clean_version/elevation/derived_for_model/global/ptha18/merged_gebco_ga250_dem_patched.tif',
    new_files, fixed=TRUE)

# Double check that they exist
new_files_exist = file.exists(new_files)
if(!all(new_files_exist)){
    print('Could not find the following files: ')
    print(new_files[which(!new_files_exist)])
    stop('Deliberate halt')
}

# Write to file
writeLines(new_files, con='swals_elevation_files_in_preference_order.txt')
