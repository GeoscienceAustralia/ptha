This folder contains sub-folders with source code:

* [shallow_water](shallow_water) contains the main codes for the shallow water equations solver
* [parallel](parallel) contains code to allow running distributed-memory parallel jobs with MPI (or coarrays)
* [raster](raster) contains code for reading from rasters via GDAL
* [util](util) contains various utilities, such as point-in-polygon, code timing, linear interpolation, etc. These do not depend heavily on the rest of the code (but everything depends on the [global_mod.f90](shallow_water/global_mod.f90) parameters).

Most modules should contain a subroutine 'test_MODULE_NAME_HERE' (where MODULE_NAME_HERE is the name of the module). This subroutine takes no arguments, but runs the modules unit tests.

In addition, there are some fragments of makefiles in the current directory, which are included in application makefiles. The core files are:
* [src_make_commands](src_make_commands) is a Makefile fragment used to compile most SWALS modules
* [src_standard_compiler_var](src_standard_compiler_var) defines key compiler variables. For better portability between different machines, this simply includes another relevant file. For example when compiling on regular Linux systems with MPI, it can include [src_standard_compiler_var_gfortran_mpi](src_standard_compiler_var_gfortran_mpi) when compiling on regular linux systems with MPI. When compiling on the Gadi supercompter on NCI, it can instead include [src_standard_compiler_var_NCI_gadi_ifort](src_standard_compiler_var_NCI_gadi_ifort).
* [test_run_commands](test_run_commands) defines key commands for running the validation tests in parallel. This lets us control the number of openmp threads, MPI ranks, etc. For better portability between different machines, this simply includes another relevant file. For example [test_run_commands_basic](test_run_commands_basic) for my home machine as of 2021, or [test_run_commands_gadi_halfnode](test_run_commands_gadi_halfnode) to run on half a node (24 cores) of NCIs Gadi supercomputer.
