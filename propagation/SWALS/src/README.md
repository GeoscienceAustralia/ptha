This folder contains sub-folders with source code:

* 'util' - various utilities, such as point-in-polygon, code timing, linear
interpolation, etc. These do not depend heavily on the rest of the code (but
everything depends on the 'global_mod.f90' parameters)
* 'raster' - interface to GDAL for reading from rasters
* 'shallow_water' - main codes for the shallow water equations solver
* 'parallel' - codes to allow running in parallel with coarrays and/or MPI

Most modules should contain a subroutine 'test_MODULE_NAME_HERE' (where MODULE_NAME_HERE is
the name of the module). This subroutine takes no arguments, and it runs the module's unit tests.

In addition, there are some fragments of makefiles in the current directory, which
are included in application specific makefiles.
