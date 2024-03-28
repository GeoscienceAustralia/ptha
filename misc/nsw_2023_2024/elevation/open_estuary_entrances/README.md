# Create polygons/values defining estuary entrances.

Run with
```
Rscript make_estuary_entrance_polygons_values_file.R
Rscript make_file_paths_absolute.R
```

The script [make_estuary_entrance_polygons_values_file.R](make_estuary_entrance_polygons_values_file.R) creates polygon,value pairs, writing outputs to a csv file with relative file paths.

The script [make_file_paths_absolute.R](make_file_paths_absolute.R) will make a new file with absolute file paths. This might need to be rerun when moving to another system.

The SWALS model will ensure the elevation inside these polygons does not exceed the specified value.
