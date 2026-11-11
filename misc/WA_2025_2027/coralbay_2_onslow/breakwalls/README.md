This directory contains a script to create x/y/z data representing the tops of breakwalls.

## Generate the breakwall data
First, generate breakwall data as you please and save it as a csv file in any subdirectory here. The csv file should have the three columns for the x/y/z (or similarly lon/lat/z) coordinates. The csv file should have a header, although it's ignored when read into SWALS.

This project makes two types of breakwalls:
1. [buffer_height/make_breakwalls.R](./buffer_height/make_breakwalls.R) generates breakwall csv files by taking a the maximum elevation in a neighbourhood around each point along a line. Edit the bottom of this R script to specify which elevation rasters to use, the spacing of points along the line, the buffer (neighbourhood) size, and optionally add a fixed offset to z.
2. [fixed_height/make_breakwalls.R](./fixed_height/make_breakwalls.R) generates a breakwall csv using a fixed elevation z along each line. Edit the bottom of this file to specify the elevation of each line.

## Generate a list of breakwall files
Then, generate a list of breakwall files to be used in the model. This is done by running the R script [make_breakwalls_list.R](make_breakwalls_list.R). This script reads the breakwall data and writes a list of breakwall files to be used in the model. It reads any csv file in any subdirectory recursively and writes a list of these files to `swals_breakwalls_list.txt`.
