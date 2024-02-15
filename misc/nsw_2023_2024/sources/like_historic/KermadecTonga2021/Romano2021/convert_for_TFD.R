#
# Convert the fault_geometry.dat (containing quadrilateral sources) into a set
# of triangles, with input format suitable for the TFD.x triangular dislocation
# code. This allows us to compute the vertical deformation from the faults
#

fault_geo = read.table('fault_geometry.dat', header=TRUE)

# For each quad, make a triangle by dropping the 4th coordinate
fault_geo_t1 = fault_geo[,-c(10, 11, 12)]
# As above, but drop the second coordinate, so the triangle covers the 'other half' of the quad
fault_geo_t2 = fault_geo[,-c(4, 5, 6)]

# Combine them
names(fault_geo_t2) = names(fault_geo_t1)
fault_geo_t = rbind(fault_geo_t1, fault_geo_t2)

## The multiseg.dat format looks like this:
# LON1	LAT1	DEPTH1(km)	LON2	LAT2	DEPTH2(km)	LON3	LAT3	DEPTH3(km)	RAKE	SLIP(m)
# 330.0000001 15.5000001 4.0 330.0000001 14.5000001 24.0 330.0000001 15.5000001 24.0 90.0 5.0

# Drop the columns we don't use.
fault_geo_t = fault_geo_t[,-c(12)]
# Make the names match the multiseg.dat format
names(fault_geo_t) = c('LON1', 'LAT1', 'DEPTH1(km)', 'LON2', 'LAT2', 'DEPTH2(km)', 'LON3', 'LAT3', 'DEPTH3(km)', 'RAKE', 'SLIP(m)')

# Subtract the min depth of any of the fault points -- assuming this represents the seafloor.
# (not done here as currently depth-->0 can lead to numerical issues in TFD).
#fault_geo_t[c('DEPTH1(km)', 'DEPTH2(km)', 'DEPTH3(km)')] = fault_geo_t[c('DEPTH1(km)', 'DEPTH2(km)', 'DEPTH3(km)')] -
#                                                       min(fault_geo_t[c('DEPTH1(km)', 'DEPTH2(km)', 'DEPTH3(km)')])

write.table(fault_geo_t, file='multiseg.dat', sep=" ", row.names=FALSE, quote=FALSE)

# Call the TFD program
system('./TFD_src_v5.0/TFD.x multiseg_derived_from_fault_geometry.dat')

# Read the deformations
xyz = scan('UZ.xyz')
dim(xyz) = c(3, length(xyz)/3)
unique_x = unique(xyz[1,])
unique_y = unique(xyz[2,])
uz = -xyz[3,]
dim(uz) = c(length(unique_x), length(unique_y))
#uz = uz[,dim(uz)[2]:1] # Bottom-up rather than top-down
dx = diff(range(unique_x))/(dim(uz)[1] - 1)
dy = diff(range(unique_y))/(dim(uz)[2] - 1)
minx = min(unique_x)-dx/2
maxx = max(unique_x)-dx/2
miny = min(unique_y)-dy/2
maxy = max(unique_y)-dy/2

# Pack into raster and save
library(terra)
outrast = rast(nrows=dim(uz)[1], ncols=dim(uz)[2], xmin=minx, xmax=maxx, ymin=miny,ymax=maxy, crs='EPSG:4326')
outrast[[1]] = as.numeric(uz)

writeRaster(outrast, 'Kermadec2021_Romano.tif', 
            gdal=c('COMPRESS=DEFLATE'), overwrite=TRUE)
