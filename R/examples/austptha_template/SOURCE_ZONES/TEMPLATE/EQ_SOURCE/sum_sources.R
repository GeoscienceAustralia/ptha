library(rasterVis)
rasts = lapply(Sys.glob('Unit_source_data/*/*.tif'), raster)
r1 = rasts[[1]]
for(i in 2:length(rasts)) r1 = r1 + rasts[[i]]

writeRaster(r1, 'sum_unit_sources.tif', options=c('COMPRESS=DEFLATE'), overwrite=TRUE)
