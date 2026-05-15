library(terra)
output_filename = 'locations.csv'

x = vect('quick_gauges_NWWA')
output_coords = data.frame(lon = crds(x)[,1], lat=crds(x)[,2], gaugeID= seq(1, length(x)))

write.csv(output_coords, output_filename, row.names=FALSE)
