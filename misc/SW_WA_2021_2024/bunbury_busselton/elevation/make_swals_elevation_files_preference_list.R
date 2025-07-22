#
# List out swals files
#

# Global DEM from PTHA18
PTHA18_DEM = '/g/data/w85/tsunami/MODELS/AustPTHA_c/DATA/ELEV/merged_dem/merged_gebco_ga250_dem_patched.tif'

# GA250m -- no we are covered by the transition DEM.
# GA250_DEM = 

# WA transition DEM
Greater_Perth_Transition_DEM = '../../../DATA/ELEV/WA/Transition_DEM_near_Perth_Oct2021/WA_smooth_between_GA250_and_MultibeamNearshore.tif'

# NW WA DEMs
NW_WA_DEMs = paste0('/g/data/w85/tsunami/DATA/ELEVATION/WA/AusSeabed_NW_Shelf/',
    c('North_West_Shelf_Satellite_Derived_Bathymetry_2020_10m_MSL_cog_WGS84.tif',
      'North_West_Shelf_DEM_v2_Bathymetry_2020_30m_MSL_cog_WGS84.tif'))

# Coastal tiles
Coastal_merge_tiles = Sys.glob('../../../DATA/ELEV/WA/SWWA_nearshore_tifs_2021/merged_rasters/*.tif')

# Patch bunbury DEM -- superceeded.
#bunbury_dem = '../../../DATA/ELEV/WA/WA_Bunbury_estuary_merged_tile13/tile13_v3_cl2.tif'

# Patch bunbury DEM -- this patches over the previous version.
bunbury_dem = '../../../DATA/ELEV/WA/WA_Bunbury_revised_November2022/Bunbury_patch_tile.tif'
# Zoom around floodgate
bunbury_floodgate_dem = '../../../DATA/ELEV/WA/WA_Bunbury_revised_November2022/Bunbury_floodgate_channel_no_bridge_piers.tif'
# Extra highres around floodgate -- this is a trick to locally circumvent SWALS's bilinear interpolation
bunbury_floodgate_dem_extrahighres = '../../../DATA/ELEV/WA/WA_Bunbury_revised_November2022/Bunbury_floodgate_only_extreme_resolution.tif'
# 2022 Onshore Busselton DEM
#busselton_1m_2022 = '/media/gareth/Windows7_OS/Users/gareth/Documents/work/AustPTHA/DATA/ELEV/WA/Busselton_files/20220330_Aerometrix_Lidar/DEM_merged_1m/Busselton_merged_20220330_Aerometrix_Lidar_1m_WGS84.tif'
busselton_1m_2022 = '../../../DATA/ELEV/WA/Busselton_files/20220330_Aerometrix_Lidar/DEM_merged_1m/Busselton_merged_20220330_Aerometrix_Lidar_1m_WGS84.tif'

# Patch DEM in the Vasse estuary, Busselton
#busselton_vasse_patch = '/media/gareth/Windows7_OS/Users/gareth/Documents/work/AustPTHA/DATA/ELEV/WA/Busselton_files/GIS_and_processing/Vasse_survey.tif'
busselton_vasse_patch = '../../../DATA/ELEV/WA/Busselton_files/GIS_and_processing/Vasse_survey.tif'

# DEM developed by Shane Martin for Busselton -- beware this isn't good in 'onshore' water-areas (channels, estuaries)
busselton_1m_SM = '/g/data/w85/tsunami/DATA/ELEVATION/WA/GA_modelling_before_2015/Busselton_StormSurge_ShaneMartin/Busselton_1m_full/Mosaic1mDEM163_WGS84_removeEdgeZeros.tif'

# Improve the bed in the Vasse Diversion drain [update for 2024]
Vasse_diversion_drain_bed_patch = '/g/data/w85/tsunami/MODELS/inundation/DATA/ELEV/WA/Vasse_diversion_drain_bed_patch/Vasse_drain_patch.tif'

# Subset of a patch-dem developed by Kaya Wilson, by patching over the coastal_merge_tiles with better data in Port Geographe. 
# This subset is focussed around Port Geographe, and represents recent changes to the harbour.
portgeographe_patch = '../../../DATA/ELEV/WA/WA_Busselton_PortGeographe_merged_tile15/WA_Busselton_PortGeographe_merged_tile15/patch_near_portgeographe.tif'

# A patch-dem for the mouth of one entrance to the Peel Harvey estuary (with less bank roughness than the coastal merge tiles)
peel_entrance_patch_dem = '../../../DATA/ELEV/WA/PeelHarvey_patch_dem/Peel_Harvey_Estuary_patch_fill.tif'

# Some offshore islands
Cocos_Keeling_tiles = paste0(
    '/g/data/w85/tsunami/DATA/ELEVATION/Cocos_Keeling_and_Christmas_Islands/AusSeabed_Christmas_Cocos_Islands/wgs84/',
    c('ausseabed_Cocos_Island_Bathymetry_2010_10m_OV.tif',
      'ausseabed_Cocos_Island_Bathymetry_2010_50m_OV.tif',
      'ausseabed_Cocos_Island_Bathymetry_2010_100m_OV.tif')
    # DELIBERATELY AVOID ausseabed_Cocos_Island_Bathymetry_2010_250m_OV.tif -- which does not join the Global DEM or other DEMs well.
    )

Christmas_Island_tiles = paste0( 
    '/g/data/w85/tsunami/DATA/ELEVATION/Cocos_Keeling_and_Christmas_Islands/AusSeabed_Christmas_Cocos_Islands/wgs84/',
    c('ausseabed_Christmas_Island_Bathymetry_2010_50m_OV.tif',
      'ausseabed_Christmas_Island_Bathymetry_2010_100m_OV.tif',
      'ausseabed_Christmas_Island_Bathymetry_250m_2010_250m_OV.tif')
)

files_in_preference_order = c(
    # Good-quality patches                              
    peel_entrance_patch_dem, 
    bunbury_floodgate_dem_extrahighres,
    bunbury_floodgate_dem,
    bunbury_dem, 
    Vasse_diversion_drain_bed_patch, # Update for 2024
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
    Greater_Perth_Transition_DEM, 
    # Global DEM
    PTHA18_DEM)

files_in_preference_order = normalizePath(files_in_preference_order)

if(all(file.exists(files_in_preference_order))){
    print('ALL FILES EXIST')
}else{
    stop('ERROR: Not all files exist')
}

# For backward compatibility, give versions with/without a file we added in 2024
k = grep('Vasse_diversion_drain_bed_patch/Vasse_drain_patch.tif', files_in_preference_order, fixed=TRUE)
stopifnot(length(k) == 1)
cat(files_in_preference_order[-k], file='swals_elevation_files_in_preference_order.txt', sep="\n")
cat(files_in_preference_order, file='swals_elevation_files_in_preference_order_2024_update.txt', sep="\n")

