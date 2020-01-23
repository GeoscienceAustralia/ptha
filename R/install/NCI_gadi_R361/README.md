Here I include the modules that are required to build and run rptha on NCI's gadi supercomputer.

The code has been revised to use the system-wide gdal install

Install like this:

    install.packages(c('sp', 'rgdal', 'rgeos', 'FNN', 'raster', 'minpack.lm', 'geometry', 'geosphere', 'rgl', 'ncdf4', 'testthat', 'devtools', 'roxygen2'))
