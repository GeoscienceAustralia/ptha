# Add some gridded hazard points
library(raster)
r1 = raster('../ELEV/merged_dem/merged_gebco_ga250_dem_patched.tif')

#
# Get other hazard points -- note the files have latitude before longitude,
# which we correct
#
haz_GA250_20 = read.table('OUTPUTS/OUTPUT_GA250_20m/haz_pts_W-180.txt', 
    skip=1, header=FALSE)
haz_GA250_20[,1:2] = haz_GA250_20[,2:1] 
haz_GA250_100 = read.table('OUTPUTS/OUTPUT_GA250_100m/haz_pts_W-180.txt', 
    skip=1, header=FALSE)
haz_GA250_100[,1:2] = haz_GA250_100[,2:1] 
haz_GA250_1000 = read.table('OUTPUTS/OUTPUT_GA250_1000m/haz_pts_W-180.txt', 
    skip=1, header=FALSE)
haz_GA250_1000[,1:2] = haz_GA250_1000[,2:1] 
haz_GEBCO2014_100 = read.table('OUTPUTS/OUTPUT_GEBCO2014_100m/haz_pts_W-180.txt', 
    skip=1, header=FALSE)
haz_GEBCO2014_100[,1:2] = haz_GEBCO2014_100[,2:1] 

# Points in GA250
existing_hp_GA250 = rbind(haz_GA250_20, haz_GA250_100, haz_GA250_1000)

# The above misses heard/macdonald, so let's add them back in
hp_keep = (haz_GEBCO2014_100[,1] > 70 & haz_GEBCO2014_100[,1] < 75 & 
    haz_GEBCO2014_100[,2] > -55 & haz_GEBCO2014_100[,2] < -50)

existing_hp = rbind(existing_hp_GA250, haz_GEBCO2014_100[hp_keep,])

# We have too many points around the other 'GA250' regions, so lets remove them
library(rgdal)
removal_polygon = readOGR(dsn='INPUTS/GRID_CLIP_LAYER/grid_clip_layer.shp', 
    layer='grid_clip_layer')
removal_polygon = removal_polygon@polygons[[1]]@Polygons[[1]]@coords
hp_remove = (point.in.polygon(existing_hp[,1], existing_hp[,2], removal_polygon[,1], removal_polygon[,2]) == 1)
nearby_hp = existing_hp[-which(hp_remove),] 

#
# Make a regular lonlat grid
#

grid_dlonlat = 1/3
g1 = expand.grid(seq(-40, 320, by=grid_dlonlat), seq(-65, 72, by=grid_dlonlat))
g1_elev = extract(r1, g1)
g1_keep = which(g1_elev < -50 & (!is.na(g1_elev)))
g1 = g1[g1_keep,]

#
# Find distance between g1 and nearby_hp, in lon/lat units
#

library(FNN)
hp_nearest_g1 = get.knnx(as.matrix(nearby_hp[,1:2]), as.matrix(g1[,1:2]), k=1)
nearest_lon = nearby_hp[hp_nearest_g1$nn.index,1]
nearest_lat = nearby_hp[hp_nearest_g1$nn.index,2]

# Remove points with lon-lat distance > 1 degree
dlonlat = 1
g1_keep = which( (abs(g1[,1] - nearest_lon) < dlonlat) & abs(g1[,2] - nearest_lat) < dlonlat)

g1 = g1[g1_keep,]
names(g1) = c('lon', 'lat')
g1 = cbind(g1, data.frame(id = 1:length(g1[,1])))

dir.create('OUTPUTS/GRIDDED_POINTS', showWarnings=FALSE)
write.csv(g1, file='OUTPUTS/GRIDDED_POINTS/gridded_points.csv', row.names=FALSE)
