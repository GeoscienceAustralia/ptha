#' Create a raster of friction values
library(terra)

working_crs <- "epsg:4326"

#  No Mangroves
lc_rast_local = rast('data/fuels_clipped.tif')
lc_rast <- project(lc_rast_local, working_crs)
lc_friction = read.csv('data/fuels_clipped_friction_no_mangroves.csv')

# Create a raster of friction values
friction_from_lc_rast <- subst(
    lc_rast,
    from=lc_friction$PSU_ID,
    to=lc_friction$manning_n,
    others=NA,
    file='manning_n_from_fuels_no_mangroves.tif',
    overwrite=TRUE
)

# cat all the friction raster normalised filenames to a text file
filename_for_swals <- 'swals_manning_n_files_in_preference_order_no_mangroves.txt'
writeLines(
    c(
        normalizePath('manning_n_from_fuels_no_mangroves.tif')
    ),
    filename_for_swals
)