#
# Replace relative file paths with absolute file paths in the estuary_entrance_polygons_values file
#
x = read.table('estuary_entrance_polygons_values_relative_file_paths.csv', header=FALSE, sep=",")
y = x
y[,1] = normalizePath(y[,1])
write.table(y, 'estuary_entrance_polygons_values.csv', row.names=FALSE, col.names=FALSE, sep=",", quote=FALSE)
