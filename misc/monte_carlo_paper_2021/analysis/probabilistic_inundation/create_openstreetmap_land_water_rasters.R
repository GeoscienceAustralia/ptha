#
# For plotting it is helpful to have rasters that show "inside-vs-outside" of the openstreetmap polygon
#

library(rptha)

tonga_coast_polygon_file = '../../elevation/Tonga_coast/Tonga_coast_nearlon180_polygon.shp'

template_rasters = Sys.glob('initial_elevation_grids/*.tif')

outdir = 'openstreetmap_land_water_grids/'
dir.create(outdir, showWarnings=FALSE)

# Make rasters that are everywhere (-1). Later we will burn "+1" into land
# regions (contained in the openstreetmap poly)
for(i in 1:length(template_rasters)){
    r1 = raster(template_rasters[i])
    r1 = r1 * 0 - 1
    writeRaster(r1, file=paste0(outdir, basename(template_rasters[i])), options=c('COMPRESS=DEFLATE'), overwrite=TRUE)
}

land_water_rasters = Sys.glob(paste0(outdir, '*.tif'))
for(i in 1:length(land_water_rasters)){
    run_command = paste0('gdal_rasterize -burn 2.0 -add -l Tonga_coast_nearlon180_polygon ', tonga_coast_polygon_file, ' ', land_water_rasters[i])
    system(run_command)
}

