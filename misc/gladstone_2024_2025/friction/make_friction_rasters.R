#' Create a raster of friction values
library(terra)

working_crs <- "epsg:4326"

## 1 
# read the land cover rast
lc_rast_local = rast('data/fuels_clipped.tif')
lc_rast <- project(lc_rast_local, working_crs)
lc_friction = read.csv('data/fuels_clipped_friction.csv')

# Create a raster of friction values
friction_from_lc_rast <- subst(
    lc_rast,
    from=lc_friction$PSU_ID,
    to=lc_friction$manning_n,
    others=NA,
    file='manning_n_from_fuels.tif',
    overwrite=TRUE
)


## 2
# Read mangrove percent cover
mangrove_friction <- read.csv('data/DEA_Mangroves_friction.csv')
dea_rast_local <- rast('data/DEA_Mangroves_Landsat_clip.tiff')
dea_rast <- project(dea_rast_local, working_crs)

# Create a raster of friction values
friction_from_dea_rast <- subst(
    dea_rast,
    from=mangrove_friction$canopy_cover_class,
    to=mangrove_friction$manning_n,
    others=NA,
    file='manning_n_from_DEA_Mangroves.tif',
    overwrite=TRUE
)

# cat all the friction raster normalised filenames to a text file
filename_for_swals <- 'swals_manning_n_files_in_preference_order.txt'
writeLines(
    c(
        normalizePath('manning_n_from_DEA_Mangroves.tif'),
        normalizePath('manning_n_from_fuels.tif')
    ),
    filename_for_swals
)