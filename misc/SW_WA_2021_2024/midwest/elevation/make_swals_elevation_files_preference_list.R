#
# List out raster files for SWALS
#

path_root <- "/g/data/w85/tsunami"

path_model_data <- paste0(path_root, "/MODELS/inundation/DATA/ELEV/WA/")
path_elev <- paste0(path_root, "/DATA/ELEVATION/")
path_data_wa <- paste0(path_root, "/DATA/ELEVATION/WA/")

# Global DEM from PTHA18
PTHA18_DEM <- paste0(path_root, "/MODELS/AustPTHA_c/DATA/ELEV/merged_dem/merged_gebco_ga250_dem_patched.tif")

# GA250m -- no we are covered by the transition DEM.
# GA250_DEM <- 

# Australia Bathymetry and Topography 2023
aust_bath_topo <- paste0(
    path_elev,
    "Australia/Australian_Bathymetry_and_Topography_2023_250m/",
    "Australian_Bathymetry_and_Topography_2023_250m_MSL_cog.tif"
    )

# WA transition DEM
Greater_Perth_Transition_DEM <- file.path(path_model_data, "Transition_DEM_near_Perth_Oct2021/WA_smooth_between_GA250_and_MultibeamNearshore.tif")

# NW WA DEMs
nw_wa_dem_files <- c("North_West_Shelf_Satellite_Derived_Bathymetry_2020_10m_MSL_cog_WGS84.tif",
      "North_West_Shelf_DEM_v2_Bathymetry_2020_30m_MSL_cog_WGS84.tif")
NW_WA_DEMs <- paste0(file.path(path_data_wa, "AusSeabed_NW_Shelf/"), nw_wa_dem_files)

# Coastal tiles
Coastal_merge_tiles <- Sys.glob(file.path(path_model_data, "SWWA_nearshore_tifs_2021/merged_rasters/*.tif"))

# Bunker Bay patch DEM
bunker_bay_patch_DEM <- file.path(path_data_wa, "Bunker_Bay/Bunker_Bay_patch_DEM.tif")

# Simple patch at Henderson (constant elevation)
henderson_patch_DEM <- file.path(path_model_data, "Henderson_patch/Henderson_patch_rast.tif")

# Patches at Rottnest
rottnest_patch_DEM <- file.path(path_data_wa, "Landgate_DEM_Rottnest_FremantlePort_for_GA_DO_NOT_DISTRIBUTE/rottnest_clip/Rottnest_Island_Feb_2023_1m_DEM_WGS84_clipped.tif")
rottnest_beach_patch_DEM <- file.path(path_data_wa, "Landgate_DEM_Rottnest_FremantlePort_for_GA_DO_NOT_DISTRIBUTE/rottnest_clip/Rottnest_beach_patch_raster.tif")
# Fremantle LiDAR which includes the reclaimed land at North Fremantle (missed by most other products)
fremantle_patch_DEM <- file.path(path_data_wa, "Landgate_DEM_Rottnest_FremantlePort_for_GA_DO_NOT_DISTRIBUTE/fremantle_clip/Fremantle_Aug2021_LiDAR_1m_DEM_WGS84_clipped.tif")


# A large-scale DEM (derived from Photogrammetry?). To prevent artefacts this was limited onshore in a crude way which is OK for our application.
kalbarri_to_isrealite_bay_photogram <- file.path(path_data_wa, "WA_Kalbarri_to_Isrealite_bay/clip_onshore/Kalbarri_to_Isrealite_bay_WGS84_clipped_onshore.vrt")

# Patch bunbury DEM -- superceeded.
#bunbury_dem <- file.path(path_model_data, "WA_Bunbury_estuary_merged_tile13/tile13_v3_cl2.tif")

# Patch bunbury DEM -- this patches over the previous version.
bunbury_dem <- file.path(path_model_data, "WA_Bunbury_revised_November2022/Bunbury_patch_tile.tif")
# Zoom around floodgate
bunbury_floodgate_dem <- file.path(path_model_data, "WA_Bunbury_revised_November2022/Bunbury_floodgate_channel_no_bridge_piers.tif")
# Extra highres around floodgate -- this is a trick to locally circumvent SWALS"s bilinear interpolation
bunbury_floodgate_dem_extrahighres <- file.path(path_model_data, "WA_Bunbury_revised_November2022/Bunbury_floodgate_only_extreme_resolution.tif")
# 2022 Onshore Busselton DEM
busselton_1m_2022 <- file.path(path_model_data, "Busselton_files/20220330_Aerometrix_Lidar/DEM_merged_1m/Busselton_merged_20220330_Aerometrix_Lidar_1m_WGS84.tif")

# Patch DEM in the Vasse estuary, Busselton
#busselton_vasse_patch <- "/media/gareth/Windows7_OS/Users/gareth/Documents/work/AustPTHA/DATA/ELEV/WA/Busselton_files/GIS_and_processing/Vasse_survey.tif"
busselton_vasse_patch <- file.path(path_model_data, "Busselton_files/GIS_and_processing/Vasse_survey.tif")

# DEM developed by Shane Martin for Busselton -- beware this isn"t good in "onshore" water-areas (channels, estuaries)
busselton_1m_SM <- file.path(path_data_wa, "GA_modelling_before_2015/Busselton_StormSurge_ShaneMartin/Busselton_1m_full/Mosaic1mDEM163_WGS84_removeEdgeZeros.tif")

# Subset of a patch-dem developed by Kaya Wilson, by patching over the coastal_merge_tiles with better data in Port Geographe. 
# This subset is focussed around Port Geographe, and represents recent changes to the harbour.
portgeographe_patch <- file.path(path_model_data, "WA_Busselton_PortGeographe_merged_tile15/WA_Busselton_PortGeographe_merged_tile15/patch_near_portgeographe.tif")

# A patch-dem for the mouth of one entrance to the Peel Harvey estuary (with less bank roughness than the coastal merge tiles)
peel_patch_dem <- file.path(path_model_data, "PeelHarvey_patch_dem/Peel_Harvey_Estuary_patch_fill.tif")

# Some offshore islands
Cocos_Keeling_tiles <- paste0(
    path_elev,
    "/Cocos_Keeling_and_Christmas_Islands/AusSeabed_Christmas_Cocos_Islands/wgs84/",
    c("ausseabed_Cocos_Island_Bathymetry_2010_10m_OV.tif",
      "ausseabed_Cocos_Island_Bathymetry_2010_50m_OV.tif",
      "ausseabed_Cocos_Island_Bathymetry_2010_100m_OV.tif")
    # DELIBERATELY AVOID ausseabed_Cocos_Island_Bathymetry_2010_250m_OV.tif -- which does not join the Global DEM or other DEMs well.
    )

Christmas_Island_tiles <- paste0( 
    path_elev,
    "/Cocos_Keeling_and_Christmas_Islands/AusSeabed_Christmas_Cocos_Islands/wgs84/",
    c("ausseabed_Christmas_Island_Bathymetry_2010_50m_OV.tif",
      "ausseabed_Christmas_Island_Bathymetry_2010_100m_OV.tif",
      "ausseabed_Christmas_Island_Bathymetry_250m_2010_250m_OV.tif")
)

# Midwest patches: rivers, estuaries and offshore areas burned into the midwest, WA
midwest_patches <- paste0(path_model_data,
    "midwest_patches/rasters/",
    c("chapman_river.tif",
     "geraldton_sheds.tif",
     "greenough_river.tif",
     "irwin_river.tif",
     "moore_river.tif",
     "off_lancelin.tif",
     "off_cocos_islands.vrt")
    )


files_in_preference_order <- c(
    # Good-quality patches
    midwest_patches,
    bunker_bay_patch_DEM, 
    henderson_patch_DEM,
    rottnest_patch_DEM,
    rottnest_beach_patch_DEM,
    fremantle_patch_DEM,
    peel_patch_dem, 
    bunbury_floodgate_dem_extrahighres,
    bunbury_floodgate_dem,
    bunbury_dem, 
    busselton_1m_2022,
    busselton_vasse_patch,
    portgeographe_patch,
    busselton_1m_SM,
    # Mainstream good-quality data
    Coastal_merge_tiles, 
    Cocos_Keeling_tiles, 
    Christmas_Island_tiles, 
    # Transition DEMs and regional products
    NW_WA_DEMs,
    kalbarri_to_isrealite_bay_photogram,
    aust_bath_topo,
    Greater_Perth_Transition_DEM,
    # Global DEM
    PTHA18_DEM)


files_in_preference_order <- normalizePath(files_in_preference_order)

file_existence <- file.exists(files_in_preference_order)
if (all(file_existence)) {
    print('ALL FILES EXIST')
} else {
    missing_files <- files_in_preference_order[!file_existence]
    stop(paste('ERROR: The following files do not exist:\n', paste(missing_files, collapse = '\n')))
}

cat(files_in_preference_order, file = "./swals_elevation_files_in_preference_order.txt", sep = "\n")
