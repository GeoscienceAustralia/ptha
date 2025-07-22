library(terra)
solomon07 = rast('../../../sources/like_historic/solomon2007/solomon2007-batch2_variable_area_uniform_slip_11244_count_1.tif')
output = solomon07 * (1/1000)

writeRaster(output, file='small_scenario_solomon2007_divided_by_1000.tif', gdal=c('COMPRESS=DEFLATE'), overwrite=TRUE)
