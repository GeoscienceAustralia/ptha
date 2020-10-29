# Make a DEM that is 0 for elevation above MSL [dry], and 1 for elevation below [wet].
# The resulting file is small so easy to distribute via THREDDS.
library(raster)
x = raster('merged_gebco_ga250_dem_patched.tif')

y = (x < 0)
crs(y) = CRS('+init=epsg:4326')

writeRaster(y, file='wet_or_dry_gebco_ga250_dem_patched.tif', options=c('COMPRESS=DEFLATE'))
