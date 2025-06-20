#
# Generate a list of elevation files in order of preferance.
# Save to file.
#

# Global DEM from PTHA18
gebco250_patched = '/g/data/w85/tsunami/MODELS/AustPTHA_c/DATA/ELEV/merged_dem/merged_gebco_ga250_dem_patched.tif'

# Australia wide Bathymetry and Topography
australia250_2023 = '/g/data/w85/tsunami/DATA/ELEVATION/Australia/Australian_Bathymetry_and_Topography_2023_250m/Australian_Bathymetry_and_Topography_2023_250m_MSL_cog.tif'

# ELVIS 1m in Gladstone area
# for rockhampton_2015, it was better to cull negative values than by WOFS
rockhampton_2015 = '/g/data/w85/tsunami/MODELS/inundation/DATA/ELEV/QLD/clip_neg/Rockhampton_2015_no_water_wgs84.tif'
gladstone_by_WOFS = paste0(
    '/g/data/w85/tsunami/DATA/ELEVATION/QLD/Gladstone_ELVIS_1m_LIDAR/wgs84_version_clipped_onshore/',
    c(  
        'StLawrence201805_DEM_wgs84/StLawrence201805_DEM_wgs84.vrt',
        'Bundaberg201801_DEM_wgs84/Bundaberg201801_DEM_wgs84.vrt',
        'Bundaberg_2016_wgs84/Bundaberg_2016_wgs84.vrt',
        'Bundaberg_2011_wgs84/Bundaberg_2011_wgs84.vrt',
        'Bundaberg_2009_wgs84/Bundaberg_2009_wgs84.vrt',
        'Gladstone201801_DEM_wgs84/Gladstone201801_DEM_wgs84.vrt',
        'Livingstone_2015_wgs84/Livingstone_2015_wgs84.vrt',
        # 'Rockhampton_2015_wgs84/Rockhampton_2015_wgs84.vrt',
        'Rockhampton201305_DEM_wgs84/Rockhampton201305_DEM_wgs84.vrt',
        'Rockhampton_2009_wgs84/Rockhampton_2009_wgs84.vrt',
        'Gladstone201307_DEM_wgs84/Gladstone201307_DEM_wgs84.vrt',
        'Gladstone_2009_wgs84/Gladstone_2009_wgs84.vrt'
    )
)
gladstone_elvis = c(rockhampton_2015, gladstone_by_WOFS)

# Rasters extracted to put in high preferance. This is to patch over small areas around the foreshore with artifacts.
raster_patches <- c(
    Sys.glob("/g/data/w85/tsunami/MODELS/inundation/DATA/ELEV/QLD/gladstone_priority_elev/priority_rasters/*.tif"),
    '/g/data/w85/tsunami/MODELS/inundation/DATA/ELEV/QLD/clip_neg/permian_point_wgs84.tif'
)


# ELVIS 1m in GBR islands near Gladstone
gladstone_islands_elvis = paste0(
    '/g/data/w85/tsunami/DATA/ELEVATION/QLD/Gladstone_ELVIS_1m_LIDAR/wgs84_version_clipped_onshore/',
    c(  
        'HeronIsland_2009_wgs84/HeronIsland_2009_wgs84.vrt',
        'LadyElliotIsland_2009_wgs84/LadyElliotIsland_2009_wgs84.vrt',
        'OneTreeIsland_2009_wgs84/OneTreeIsland_2009_wgs84.vrt'
    )
)

# Gladstone Port Authority: post drege data
gladstone_drege = '/g/data/w85/tsunami/DATA/ELEVATION/QLD/Gladstone_Post_Dredge/convert_to_wgs84/GL020119_Post-Dredge_0p5m_Cube_wgs84_AHD.tif'

# Gladstone Port Corporation
port_alma_channel = '/g/data/w85/tsunami/DATA/ELEVATION/QLD/Gladstone_port_corporation/data_clean/pa_may2021.tif'
fishermans_landing = '/g/data/w85/tsunami/DATA/ELEVATION/QLD/Gladstone_port_corporation/data_clean/wbra_jfp_2019.tif'

# Gladstone MSG data
gladstone_msg = c(
    '/g/data/w85/tsunami/DATA/ELEVATION/QLD/MSG_Gladtstone_Port/Gladstone_Data_full_density/GL020120_Gladstone_Historic_data_Interpolated_Merged_2p5m_Cube_FINAL_AHD_wgs.tif',
    '/g/data/w85/tsunami/DATA/ELEVATION/QLD/MSG_Gladtstone_Port/Gladstone_Data_upto_5m_contour/MSG_Gladstone_port_10m_wgs84.tif'
)

# Gladstone rivers
gladstone_rivers = paste0('/g/data/w85/tsunami/DATA/ELEVATION/QLD/Gladstone_Regional_Council/',
    c(
        'Auckland_Creek_DEM/AucklandCk_5m_bathy_wgs84.tif',
        'Auckland_Creek_DEM/AucklandCk_Marina_5m_bathy_wgs84.tif',
        'Baffle_Creek_Bathymetric_Raster/Baffle_Creek_Bathymetric_Survey_by_centreline/Baffle_Creek_Bathymetric_Survey_by_centreline.vrt',
        'Boyne_Bathymetric_Data_wgs84.tif'
    )
)

# NSW LIDAR from ELVIS
nsw_lidar = '/g/data/w85/tsunami/DATA/ELEVATION/NSW/coastal_lidar/NSW_5m_near_coast/tifs_masked_by_WOFS/NSW_onshore_Lidar_5m_masked_by_WOFS.vrt'

# Great Barrier reef
barrier_reef = '/g/data/w85/tsunami/DATA/ELEVATION/GreatBarrierReef_30m/Great_Barrier_Reef_30m_2020/GBR_2020_merged.vrt'

# Patches for this model
flat_patches <- Sys.glob("/g/data/w85/tsunami/MODELS/inundation/DATA/ELEV/QLD/gladstone_patches/flat_rasters/*.tif")

top_flat_patches <- Sys.glob("/g/data/w85/tsunami/MODELS/inundation/DATA/ELEV/QLD/gladstone_patches/top_flat_rasters/*.tif")


files_in_preference_order = c(
    # Better than a patch
    fishermans_landing,
    # Patches
    top_flat_patches,
    raster_patches,
    flat_patches,
    # Good-quality rasters
    gladstone_drege,
    port_alma_channel,
    gladstone_rivers,
    gladstone_msg,
    # Mainstream good-quality data
    gladstone_elvis,
    gladstone_islands_elvis,
    nsw_lidar,
    # Transition DEMs and regional products
    barrier_reef,
    australia250_2023,
    # Global DEM
    gebco250_patched)

files_in_preference_order = normalizePath(files_in_preference_order)

if(!all(file.exists(files_in_preference_order))){
    stop('ERROR: Not all files exist')
}

# Create the file for swals 
cat(files_in_preference_order, file='swals_elevation_files_in_preference_order.txt', sep="\n")

# Create links to the data directories, can be useful for locating and copying raw data 
path_data = "./data_links"
dir.create(path_data, showWarnings = FALSE)
for (i in seq_along(files_in_preference_order)) {
    # if the symbolic link exists skip
    if (file.exists(paste0(path_data, "/", basename(dirname(files_in_preference_order[i]))))) {
        next
    } else{
        ln_command = paste0("ln -s ", dirname(files_in_preference_order[i]), " ", path_data)
        system(ln_command)
    }
}

current_time = format(Sys.time(), "%a %b %d %X %Y")
cat(current_time, file="swals_elevation_files_in_preference_order.log", append=FALSE)
