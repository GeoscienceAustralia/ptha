# Specify elevation rasters used by the model. 

The raster filenames are in `swals_elevation_files_in_preference_order.txt`.
This file was generated from QGIS by ordering the raster layers into a group called `Raw Elevation` and then running 
[save_qgis_layers.py](save_qgis_layers.py) from the python console within QGIS.

The folder [limit_elevation_maxima_polygons](./limit_elevation_in_polygons/) is used to limit the elevation maxima in some regions. Call it using:
```{bash}
. ../../modules_R_431.sh
Rscript create_elevation_limit_csv.R
```
