#
# Ho et al. invert the free-surface displacement in terms of a basis of
# 'lon-lat' gaussians.
#

lonlat = expand.grid(
    seq(360-78, 360-70, by=2/60),
    seq(-47, -36, by=2/60))

inverted_amp = read.table('jgrb53314-sup-0001-2018jb016996-si.txt', skip=10)
names(inverted_amp) = c('lon', 'lat', 'amp')

gaussian_fun<-function(lon, lat, base_lon, base_lat, sigma=15/60){
    # The basis functions look like this -- with 15-min standard deviation.
    exp(- ( (lon-base_lon)^2 + (lat-base_lat)^2 )/(2 * sigma^2))
}

z = lonlat[,1]*0
for(i in 1:length(inverted_amp$lon)){
    z = z + inverted_amp$amp[i] * gaussian_fun(lonlat[,1], lonlat[,2],
        base_lon=inverted_amp$lon[i], base_lat=inverted_amp$lat[i])
}

library(raster)
free_surface = rasterFromXYZ(cbind(lonlat, z), crs=CRS("+init=epsg:4326"))

writeRaster(free_surface, file='Ho_Chile1960_initial_displacement.tif', options=c('COMPRESS=DEFLATE'))
