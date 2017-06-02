# This script can make downdip lines from the source contours, given
# a desired unit source length.
#
# Typically the downdip lines would be used to define the unit source
# lateral boundaries, optionally with further manual editing

source_shapefile = 'SOURCEZONE_CONTOURS/sunda.shp'
out_shapefile = 'SOURCEZONE_DOWNDIP_LINES/sunda_downdip.shp'


desired_unit_source_length = 50

library(rptha)

source_contours = readOGR(source_shapefile, 
    layer=gsub('.shp','',basename(source_shapefile)))

ds1 = create_downdip_lines_on_source_contours_improved(
    source_contours, 
    desired_unit_source_length=desired_unit_source_length)

out_shp = downdip_lines_to_SpatialLinesDataFrame(ds1)

writeOGR(out_shp, 
    dsn=out_shapefile, 
    layer=gsub('.shp', '', basename(out_shapefile)),
    driver='ESRI Shapefile', overwrite=TRUE)

