Here I include the modules that are required to build and run rptha on NCI's gadi supercomputer.

Note that I'm using a private gdal build -- the only reason for this is that at the time of writing, NCI hadn't installed gdal on gadi - but in theory a standard module should be built, and it should work fine.

Interestingly this time, the package installs worked as per the standard instructions:

    install.packages(c('sp', 'rgdal', 'rgeos', 'FNN', 'raster', 'minpack.lm', 'geometry', 'geosphere', 'rgl', 'ncdf4', 'testthat', 'devtools', 'roxygen2'))
