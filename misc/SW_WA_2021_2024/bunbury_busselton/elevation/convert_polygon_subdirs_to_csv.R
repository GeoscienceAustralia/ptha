library(sf)

# Vasse estuary polygon
vasse_estuary = st_read('initial_stage_40cmAHD/initial_stage_40cmAHD.shp')
poly_coords = st_coordinates(vasse_estuary$geometry[1])
write.csv(data.frame(lon=poly_coords[,1], lat=poly_coords[,2]), 
    file='initial_stage_40cmAHD/initial_stage_40cmAHD.csv', row.names=FALSE)

# Bridges
bridges = st_read('bridges_to_remove/bridges_to_remove.shp')
for(i in 1:length(bridges$geometry)){
    poly_coords = st_coordinates(bridges$geometry[i])
    write.csv(data.frame(lon=poly_coords[,1], lat=poly_coords[,2]), 
        file=paste0('bridges_to_remove/bridges_to_remove_', i, '.csv'), 
        row.names=FALSE)
}
output_files = paste0('../elevation/bridges_to_remove/bridges_to_remove_', seq(1, length(bridges$geometry)), '.csv')
cat(output_files, file='force_elevation_to_zero_or_below_files.txt', sep="\n")
