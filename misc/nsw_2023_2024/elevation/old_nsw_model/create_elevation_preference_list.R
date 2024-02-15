# Various elevation datasets 

# NSW 2018 LIDAR/LADS near the coast
bathy_2018 = Sys.glob("bathy2018_patched/*.tif")

# A few other good DEMs that come second preference
intermediate_improvements = c(Sys.glob("good_patch_dems/*.tif"),
    "orig/vic/vcdem2017/Victorian-coast_Bathy_10m_EPSG4326.tif")

# Third preference -- onshore coastal lidar [note we tried to remove water areas and artefacts]
coastal_lidar = Sys.glob("coastal_lidar_masked_by_WoFS_patched/*.tif")

# Fourth preference -- interpolated point-data in the estuaries from OEH singlebeam archive. This data
# has loads of artefacts because of our naive triangulation, but in the majority of cases the problematic
# areas are covered by higher preference datasets.
estuary_patches = Sys.glob("Gap_filling_estuaries/all_patches/*.tif")

# Some other datasets that fill remaining gaps
lesser_improvements = Sys.glob("lesser_patch_dems/*.tif")

# 30m Great barrier reef DEM
gbr30m = 'orig/qld/gbr30_all/dblbnd.adf'

# Transition DEMs for NSW and Victoria
ocean_smooth_transition = c(
    "ocean_smooth_transition_east_coast/transition_DEM_151.4_154_-33.25_-28.9.tif",
    "ocean_smooth_transition/transition_DEM_149_153_-38_-33.25.tif",
    "derived_for_model/vic/smooth_between_GA250m_and_coast/Victoria_smooth_between_GA250_and_Vic10m.tif")

# Background DEM
ptha18_global_dem = 'derived_for_model/global/ptha18/merged_gebco_ga250_dem_patched.tif'


#
# Write a text file that lists them all in order of preference
# 

all_files = c(bathy_2018, intermediate_improvements, coastal_lidar, estuary_patches, 
    lesser_improvements, gbr30m, ocean_smooth_transition, ptha18_global_dem)
# Give the full file path in all cases
all_files_normalised = normalizePath(all_files)

cat(all_files_normalised, file="swals_elevation_files_in_preference_order.txt", sep="\n")

