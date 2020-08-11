#
# Convert the fault_geometry.dat (containing quadrilateral sources) into a set
# of triangles, with input format suitable for the TFD.x triangular dislocation
# code. This allows us to compute the vertical deformation from the faults
#

fault_geo = read.table('fault_geometry.dat', header=TRUE)

# Path to the TFD.x executable in TFD_src_v4.0. GD obtained this from Fabrizio
# Romano, who re-wrote from Mead (2007) and used several other codes.
path_to_TFD.x = '/media/gareth/Windows7_OS/Users/gareth/Documents/work/AustPTHA/DATA/INVERSIONS/fabrizioRomano/Sources4Gareth/TFD_src_v4.0/TFD.x'
if(!file.exists(path_to_TFD.x)){
    stop('Could not find the TFD.x software executable at the specified path. Assuming you have compiled the code, you need to correct the path.')
}


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
fault_geo_t = fault_geo_t[,-c(10, 11, 14)]
# Make the names match the multiseg.dat format
names(fault_geo_t) = c('LON1', 'LAT1', 'DEPTH1(km)', 'LON2', 'LAT2', 'DEPTH2(km)', 'LON3', 'LAT3', 'DEPTH3(km)', 'RAKE', 'SLIP(m)')

# Subtract the min depth of any of the fault points -- assuming this represents the seafloor.
# (not done here as currently depth-->0 can lead to numerical issues in TFD).
#fault_geo_t[c('DEPTH1(km)', 'DEPTH2(km)', 'DEPTH3(km)')] = fault_geo_t[c('DEPTH1(km)', 'DEPTH2(km)', 'DEPTH3(km)')] -
#                                                       min(fault_geo_t[c('DEPTH1(km)', 'DEPTH2(km)', 'DEPTH3(km)')])

write.table(fault_geo_t, file='multiseg_derived_from_fault_geometry.dat', sep=" ", row.names=FALSE, quote=FALSE)

system(paste0(path_to_TFD.x, ' multiseg_derived_from_fault_geometry.dat'))

# Read the deformations
ux = matrix(scan('UX.xyz'), ncol=3, byrow=TRUE)
uy = matrix(scan('UY.xyz'), ncol=3, byrow=TRUE)
uz = matrix(scan('UZ.xyz'), ncol=3, byrow=TRUE)

plot(c(fault_geo_t$LON1, fault_geo_t$LON2, fault_geo_t$LON3), 
     c(fault_geo_t$LAT1, fault_geo_t$LAT2, fault_geo_t$LAT3), 
     asp=1/cos(median(fault_geo_t$LAT1)/180*pi), pch=19)
 
for(i in 1:nrow(fault_geo_t)){
   polygon(fault_geo_t[i,c(1,4,7, 1)], fault_geo_t[i,c(2,5,8,2)], 
           col=c('white', rev(rainbow(30)))[ceiling(fault_geo_t[i,11])], border='grey')
}

library(raster)
uplift = rasterFromXYZ(uz, digits=2)
plot(uplift)
crs(uplift) = CRS("+init=epsg:4326")
writeRaster(uplift, 'Illapel_2015_Romano.tif', options=c('COMPRESS=DEFLATE'), overwrite=TRUE)
