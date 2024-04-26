# Make an updated representation of the Vasse Diversion Drain walls in a particular region, to
# account for some walls built in 2020 that were not otherwise in our model or input data.
# The inputs were defined using data for the walls provided by WA water coorporation.
#
# This will serve as a new "minimum elevation" in the model, used IN ADDITION to the existing breakwalls.
#
#

library(terra)
wall_elevations = vect('distances_initial.shp')
wall_xyz = cbind(crds(wall_elevations), wall_elevations$elev)
colnames(wall_xyz) = c('lon', 'lat', 'z')
wall_xyz = as.data.frame(wall_xyz)

latrange = range(wall_xyz$lat)


# Get previously digitized points for the east bank
points_east_bank = read.csv('westBusselton3.csv')
keep = which(points_east_bank$lat > latrange[1] & points_east_bank$lat < latrange[2])
points_east_bank = points_east_bank[keep,]
points_east_bank$z = approx(wall_xyz$lat, wall_xyz$z, xout=points_east_bank$lat)$y

write.csv(points_east_bank, 'VasseDiversionDrainWallPatch_east.csv', row.names=FALSE)

points_west_bank = read.csv('westBusselton1.csv')
keep = which(points_west_bank$lat > latrange[1] & points_west_bank$lat < latrange[2])
points_west_bank = points_west_bank[keep,]
points_west_bank$z = approx(wall_xyz$lat, wall_xyz$z, xout=points_west_bank$lat)$y

write.csv(points_west_bank, 'VasseDiversionDrainWallPatch_west.csv', row.names=FALSE)
