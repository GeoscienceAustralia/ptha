#
# Define estuary entrances
# The code here does not treat 'holes' in estuaries
#
library(terra)

# INPUTS
x = vect('open_estuaries_2023_11_14/open_estuaries.shp')
open_entrance_value = 0.2
output_dir = 'open_entrance_polygons'
# END INPUTS

dir.create(output_dir, showWarnings=FALSE)

x_geo = geom(x, df=TRUE)
if(any(x_geo$hole != 0)){
    stop("Found a polygon with an internal hole -- this isn't yet supported")
}

# Separate out each polygon
x_geom_parts = paste0(as.character(x_geo$geom), '_', as.character(x_geo$part))
unique_geom_parts = unique(x_geom_parts)

# Write each polygon to a separate csv
output_files = rep(NA, length(unique_geom_parts))
output_files_values = rep(NA, length(unique_geom_parts))
for(i in 1:length(unique_geom_parts)){
    k = which(x_geom_parts == unique_geom_parts[i])

    output_file = paste0(output_dir, '/estuary_entrance_', unique_geom_parts[i], '.csv')

    local_geo = data.frame(lon=x_geo$x[k], lat=x_geo$y[k])

    write.csv(local_geo, output_file, row.names=FALSE)

    output_files[i] = output_file
    output_files_values[i] = open_entrance_value

}

write.table(data.frame(files=output_files, values=output_files_values), 
    file='estuary_entrance_polygons_values_relative_file_paths.csv', 
    row.names=FALSE, col.names=FALSE, sep=",", quote=FALSE)
