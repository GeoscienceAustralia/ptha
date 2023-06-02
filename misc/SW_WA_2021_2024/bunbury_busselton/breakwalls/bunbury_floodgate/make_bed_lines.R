#
# hack to ensure we get the 'gap' in the floodwall
#

library(sf)
MAX_BED_ELEV = -2.2
pts = st_read('TEMP_POINTS_BED.shp')

output = cbind(st_coordinates(pts), data.frame(elev=rep(MAX_BED_ELEV, length(pts)) ))
names(output) = c('lon','lat', 'z')

write.csv(output, file='bunbury_floodgate_bed_enforcement.csv', row.names=FALSE, quote=FALSE)
