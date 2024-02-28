This directory contains a script to create x/y/z data representing the bottoms of channels.

These can be 'burned' into the elevation grid in SWALS to ensure that flow paths
are represented, even if sometimes the model is too coarse to reliably
capture the the minima.

The code is very similar to that in ../breakwalls. Although no inverts are used in this model.

First, create the inverts as line shapefiles in subdirectories. The function `make_lines_with_min_buffer` inside [make_inverts_perth.R](make_inverts_perth.R) will then read these and create a list of inverts to be used in the model. The script also contains an example that find any shapefile matching `*/*.shp`.

Then generate the full file paths for these using `Rscript make_inverts_list_have_fullpath.R`.
