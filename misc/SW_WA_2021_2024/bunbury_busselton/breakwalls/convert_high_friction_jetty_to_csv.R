library(sf)
jetty = st_read('high_friction_jetty/high_friction_jetty.shp')
poly_coords = st_coordinates(jetty$geometry[1])
write.csv(data.frame(lon=poly_coords[,1], lat=poly_coords[,2]), 
    file='high_friction_jetty/high_friction_jetty.csv', row.names=FALSE)

