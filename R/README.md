rptha
-----

rptha is the main workhorse R package. To build it, go inside 'rptha', start R, and do:

    source('build_package.R')

This will make an R package file in the directory above the package, which can be installed on the command line with:

    R CMD INSTALL rptha_XXXXX.tar.gz

where the XXXX are adapted to match the file name (and on Ubuntu, you would add 'sudo' to the start).

If the above fails because you are missing packages, then try running this prior to the install to get the required packages:

    install.packages(c('sp', 'rgdal', 'rgeos', 'FNN', 'raster', 'minpack.lm', 'geometry', 'geosphere', 'rgl', 'testthat', 'devtools'))

If you are separately reading ncdf rasters, you will also need the 'ncdf4' package.

source_contours_2_unit_sources
------------------------------

This contains example code to make tsunami unit sources from source contours


make_hazard_points
------------------

Example scripts to make hazard points (i.e. offshore points where the tsunami tide-gauges are recorded)


combine_tsunami_sources
-----------------------
Example script to convert from tsunami unit sources to tsunami event offshore
wave heights (assuming the URSGA solver was used for tsunami propagation)
