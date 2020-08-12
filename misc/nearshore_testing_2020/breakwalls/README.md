# Breakwalls

This directory contains a script to create x/y/z data representing the tops of
breakwalls, as well as the result of running that script

These can be 'burned' into the elevation grid in SWALS to ensure that barriers
to flow are represented, even if sometimes the model is too coarse to reliably
capture the breakwall maxima.

## How to run it

To make the xyz data using the shapefiles, do this:

    Rscript make_breakwalls.R

The code will write an csv file inside each subdirectory (one per shapefile), and these
are referenced by the SWALs model to ensure the breakwalls are burned into the grid.
