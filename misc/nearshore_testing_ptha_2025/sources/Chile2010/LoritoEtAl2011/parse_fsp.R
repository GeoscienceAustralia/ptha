# Convert the *.fsp file into something easier to work with
# These calculations are cross-checked by examining the supplementary table in the original paper.
## Lorito, S., F. Romano, S. Atzori, X. Tong, A. Avallone, J. McCloskey, M. Cocco, E. Boschi, and A. Piatanesi (2011). Limited overlap between the seismic gap and coseismic slip of the great 2010 Chile earthquake, Nature Geoscience, 4, 173–177, doi:10.1038/ngeo1073 

FSP = readLines('s2010MAULEC02LORI.fsp')

# Much of the tablular data is on lines that do not begin with '%'
table_data_1 = read.table(text=FSP[grep('^ ', FSP)], skip=1)
names(table_data_1) = c('lat', 'lon', 'x', 'y', 'depth', 'slip')
# Note lon/lat/depth refer to the TOP-CENTRE in this file -- different to the supplementary information in the paper.

# We still need the 'strike', 'dip' and 'rake'.
d2 = FSP[grep(' DIP =', FSP)]
strike = sapply(d2, f<-function(x) as.numeric(strsplit(x, ' ')[[1]][10]), USE.NAMES=FALSE)
dip = sapply(d2, f<-function(x) as.numeric(strsplit(x, ' ')[[1]][21]), USE.NAMES=FALSE)
# Rake doesn't seem to be in the FSP file? There is an average rake (at the top), but actually
# the inversion in the paper uses 3 different rake values
rake = c(rep(110, 64), rep(105, 144-64), rep(120, 200-144))
width = rake*0 + 25
length = rake*0 + 25
slip = table_data_1$slip
# Convert depth to centroid
depth_c = table_data_1$depth + width/2 * sin(dip/180*pi)

library(geosphere)
lonlat_c = destPoint(cbind(table_data_1$lon, table_data_1$lat), 
                     b=strike+90, # Move in this bearing
                     d = width/2 * cos(dip/180*pi) * 1000, # Move this distance in m 
                     f=0) # Spherical earth (f=0)
lonlat_c[,1] = 360 + lonlat_c[,1]

output = data.frame(lon_c = lonlat_c[,1], lat_c = lonlat_c[,2], depth_c = depth_c, 
                    length=length, width=width, slip=slip, strike=strike, dip=dip, rake=rake)
write.csv(output, file='LoritoEtAl2011_source_model.csv', row.names=FALSE)
