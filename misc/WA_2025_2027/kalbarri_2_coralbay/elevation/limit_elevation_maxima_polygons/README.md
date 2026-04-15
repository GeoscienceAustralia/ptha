# Create polygons/values defining upper limits on the elevation.

The SWALS model will ensure the elevation inside these polygons does not exceed the specified value.


Run with
```
# Combine manual shapefile geometries in "limit_elevation_polys"
# with geometries made elsewhere.
# Output to limit_elevation_polys_edited.
Rscript edit_limit_elevation_polys.R

# Create polygon files that SWALS can use
Rscript make_estuary_entrance_polygons_values_file.R 

# Ensure SWALS can read the files nomatter where it is running.
Rscript make_file_paths_absolute.R 
```
