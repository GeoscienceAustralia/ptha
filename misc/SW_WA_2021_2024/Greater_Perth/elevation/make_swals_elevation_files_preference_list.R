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

# Kaya's Bunbury DEM
bunbury_dem = '../../../DATA/ELEV/WA/WA_Bunbury_estuary_merged_tile13/tile13_v3_cl2.tif'

# FIXME: Add in Kaya's Busselton DEM, which makes a guess at the estuary bathymetry (because we don't have proper data).
# /g/data/w85/tsunami/MODELS/inundation/DATA/ELEV/WA/WA_Busselton_PortGeographe_merged_tile15


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

#print('Removing Bunbury DEM for now, please add it back when ready')

files_in_preference_order = c(peel_entrance_patch_dem, 
    bunbury_dem, 
    Coastal_merge_tiles, Cocos_Keeling_tiles, Christmas_Island_tiles, 
    NW_WA_DEMs,
    Greater_Perth_Transition_DEM, PTHA18_DEM)

files_in_preference_order = normalizePath(files_in_preference_order)

if(all(file.exists(files_in_preference_order))){
    print('ALL FILES EXIST')
}else{
    stop('ERROR: Not all files exist')
}

cat(files_in_preference_order, file='swals_elevation_files_in_preference_order.txt', sep="\n")
