
in_text = readLines('dart_historical/dartmeta_public.xls')

# Lines 2 and the last line are missing. This identifies them
small_lines = which(nchar(in_text, type='bytes') < 50)
in_text = in_text[-small_lines]

# Break up into rows
in_text_split = strsplit(in_text, split='\t', fixed=TRUE, useBytes=TRUE)

# Hack it into a data.frame
in_text_mat = matrix(unlist(in_text_split), ncol=25, byrow=TRUE)
dart_dat = data.frame(in_text_mat[-1,], stringsAsFactors=FALSE)
names(dart_dat) = in_text_mat[1,]

# Lat/Long need more work to parse to numeric
lon_numerals = gsub('\\D', ' ', dart_dat$'Buoy Longitude')
lon_ew = c('E', 'W')[1+grepl('W', enc2native(dart_dat$'Buoy Longitude'))]
lon = read.table(text=lon_numerals, sep=" ")
lon = (-1 + 2*(lon_ew == 'E'))*(as.numeric(lon[,1]) + as.numeric(lon[,3])/60 + as.numeric(lon[,5])/(60*60))
lat_numerals = gsub('\\D', ' ', dart_dat$'Buoy Latitude')
lat_ns = c('N', 'S')[1+grepl('S', enc2native(dart_dat$'Buoy Latitude'))]
lat = read.table(text=lat_numerals, sep=" ")
lat = (-1 + 2*(lat_ns == 'N'))*(as.numeric(lat[,1]) + as.numeric(lat[,3])/60 + as.numeric(lat[,5])/(60*60))

dart_dat = cbind(data.frame(lon=lon, lat=lat), dart_dat)
dir.create('dart_metadata', showWarnings=FALSE)
write.csv(dart_dat, 'dart_metadata/dart_metadata.csv', row.names=FALSE)

library(rgdal)
coordinates(dart_dat) = cbind(lon, lat)
proj4string(dart_dat) = "+init=epsg:4326"

# Save as shapefile
writeOGR(dart_dat, dsn='dart_locations', layer='dart_locations', driver='ESRI Shapefile', overwrite=TRUE)

