# Elevation raster file paths

The file `swals_elevation_files_in_preference_order.txt` should contain the raster files to use in the hydrodynamic model, with absolute file paths, ordered from high to low preference.
* On NCI that file is created by running `create_elevation_preference_list_NCI.R`
* On my home machine I used `create_elevation_preference_list_home_machine.R`. 

Originally the file list was created using GIS to overlap the rasters, and identify a good combination.

To double check that the files on NCI are actually the same as those on my home machine (where the initial inspections were done) I used `compare_md5sum_NCI_vs_home.R`. 

## Open estuary entrances

The folder `open_estuary_entrances` sets up files to limit the maximum elevation in estuary entrances that have been manually defined with polygons.

