This folder contains sub-folders with source code:

* [shallow_water](shallow_water) contains the main codes for the shallow water equations solver
* [parallel](parallel) contains code to allow running distributed-memory parallel jobs with MPI (or coarrays)
* [raster](raster) contains code for reading from rasters via GDAL
* [util](util) contains various utilities, such as point-in-polygon, code timing, linear interpolation, etc. These do not depend heavily on the rest of the code (but everything depends on the [global_mod.f90](shallow_water/global_mod.f90) parameters).

Most modules should contain a subroutine 'test_MODULE_NAME_HERE' (where MODULE_NAME_HERE is the name of the module). This subroutine takes no arguments, but runs the modules unit tests.

In addition, there are some fragments of makefiles in the current directory, which are included in application makefiles.
