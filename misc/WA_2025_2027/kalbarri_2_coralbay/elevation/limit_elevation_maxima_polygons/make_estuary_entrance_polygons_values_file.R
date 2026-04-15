#
# Define polygons where the elevation cannot exceed the prescribed limit.
#
# Useful for managing elevation artefacts, burning channels, etc.
#
# The code here does not work for polygons with 'holes'.
#
library(terra)

# INPUTS
x = vect('limit_elevation_polys_edited/limit_elevation_polys_edited.shp')
output_dir = 'limit_elevation_polys_for_swals'
# END INPUTS

dir.create(output_dir, showWarnings=FALSE)

store_output_file = rep("", length(x))
store_limit = rep(NA, length(x))
for(i in 1:length(x)){

    poly_coords = crds(x[i,])
    poly_limit = x$limit[i]

    output_file = paste0(output_dir, '/polygon_', i, '_with_limit_', poly_limit, '.csv')

    store_output_file[i] = output_file
    store_limit[i] = poly_limit

    colnames(poly_coords) = c('lon', 'lat')
    write.csv(poly_coords, file=output_file, row.names=FALSE)
}

write.table(data.frame(files=store_output_file, upper_limit_elev=store_limit), 
    file='limit_elevation_polygons_values_relative_file_paths.csv', 
    row.names=FALSE, col.names=FALSE, sep=",", quote=FALSE)
