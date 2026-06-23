#
# Compute the tidal adjustment by differencing the adjusted/unadjusted rasters
#
# BEWARE: The SWALS model prevents linear/leapfrog_linear_with_nonlinear_friction
# solvers from having elevation between "ambient_sea_level -1" and "ambient_sea_level",
# and that is applied AFTER the tidal adjustment. In such areas, the difference between
# the elevation rasters will not really be the tidal adjustment. However, this should only affect
# "global grid priority domain" cells that were adjusted.
#

library(terra)

adjusted_rasters = Sys.glob('../elevation_in_model_with_tidal_adjustment/*.tif')

unadjusted_rasters = gsub('_with_tidal_', '_no_tidal_', adjusted_rasters)
stopifnot(all(file.exists(unadjusted_rasters)))

for(i in 1:length(adjusted_rasters)){

    r1 = rast(adjusted_rasters[i])
    r2 = rast(unadjusted_rasters[i])
    r3 = r2 - r1

    output_raster = gsub('elevation0_', 'tidal_adjustment_', basename(adjusted_rasters[i]))
    writeRaster(r3, output_raster, gdal=c('COMPRESS=DEFLATE'), overwrite=TRUE)

}

system('gdalbuildvrt -resolution highest all_tidal_adjustment.vrt tidal_adjust*.tif')
