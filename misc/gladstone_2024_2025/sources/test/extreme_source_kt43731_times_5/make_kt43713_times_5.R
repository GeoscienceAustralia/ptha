library(terra)
kt43731 = rast(
    '/g/data/w85/tsunami/MODELS/inundation/east_australian_coast_2021_02/sources/kt_multi_site_selection/kermadectonga2_stochastic_43731_9.4_10000.tif'
)
output = kt43731 * 5

writeRaster(output, file='extreme_scenario_kt43731_times_5.tif', gdal=c('COMPRESS=DEFLATE'), overwrite=TRUE)
