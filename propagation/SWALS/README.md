SWALS
-----

SWALS includes a number of "Shallow WAter Like Solvers". These solve the linear
and nonlinear shallow water equations (the latter with Manning friction) using
either Cartesian or Spherical Coordinates. 

Installation prerequisites
--------------------------

## Fortran compiler
Compilation requires a version of GNU-Fortran or ifort, although the test-suite
is setup for gfortran by default. The author hasn't comprehensively tested various
compiler versions, but it is likely that gfortran versions 5 or greater will
work (without coarrays), and ifort versions 19 and greater. 

The code can run in parallel with both shared memory (openmp) and distributed
memory approaches. Distributed parallel support requires that MPI is installed,
and/or that the compiler has adequate coarray support. On gfortran, coarrays
work well using the OpenCoarrays library. Because compiler support for coarrays
is generally less mature than MPI, we have also created options for running the
code with MPI (i.e replacing coarray communication with MPI). In general 
all these parallel approaches can be used together.

## C compiler
The code also makes limited use of C to interface with the GDAL library (thus
providing flexible raster-based input) and qsort. Thus a C compiler is
required, typically gcc.

## Make
The "make" program is used for compilation.

## NetCDF and GDAL
"NetCDF" and "GDAL" must also be installed. The command-line programs
"gdal-config" and "nc-config" or "nf-config" should be available, because the
makefiles use them to identify the associated libraries. The netcdf install
should include Fortran support, and this requires building netcdf with the same
Fortran compiler as you use to build SWALS.

## ncview (optional)
The program 'ncview' is useful for quickly viewing animations of the netcdf
outputs - but this is not essential.

## R (to run the validation tests)
The validation tests and associated plots make heavy use of bash shell scripts,
and the open-source language R (to post-process model outputs), both of which
need to be installed to run the validation tests. Note R is not required to
run SWALS itself. Any other language could be used similarly to process the
output netcdf files. 

Our validation test R scripts are also reliant on a number of packages. These
can be installed from inside R with the following command:
    
    # Start R and then execute this command
    install.packages(c('ncdf4', 'fields', 'raster', 'sp', 'rgdal', 'maptools', 'rasterVis'))


Compiling and testing
---------------------

The standard build scripts are based on
[src_standard_compiler_var](./src/src_standard_compiler_var) and
[src_make_commands](./src/src_make_commands). These scripts are 'included' in
application specific makefiles.  

The script [src_standard_compiler_var](./src/src_standard_compiler_var) simply
points to another compiler specific script (there are variants for both
[gfortran](./src/src_standard_compiler_var_gfortran) and
[ifort](./src/src_standard_compiler_var_ifort) and for 
[the gadi machine on NCI](./src/src_standard_compiler_var_NCI_gadi_ifort). Variables in these scripts can
be overridden by defining them in the application-specific makefile (see the
examples). This is required in many situations (e.g. to use spherical
coordinates, different compilers or compiler options, or non-standard library locations.)

See the validation test suite for examples, and documentation of compiler options below.

# Step 1: Run the unit-tests

The unit test suite is in [./tests/unit_tests](./tests/unit_tests). These are most
useful for ensuring that your install is working, and to help us avoid
accidentally breaking features as the code evolves. When trying to install
SWALS, the first thing to do is try to compile and successfully run the unit
tests. To do this, open a terminal in that directory and do

    make -B -f make_test
    ./unit_tests > outfile.log

This should take a minute, and write a lot of "PASS" statements to outfile.log. There
should not be any "FAIL" statements or other warning messages. Open up outfile.log to check.
You can count PASS or FAIL statements like this:

    # Count the passes
    grep "PASS" outfile.log | wc -w
    # Count the failures (should be 0)
    grep "FAIL" outfile.log | wc -w

# Step 2: Run the parallel unit-tests

The parallel unit-test suite is in [./tests/parallel_tests](./tests/parallel_tests). Note
this only tests the coarray parallel operations, not OpenMP. To run the tests, open
a terminal in that directory and do:

    # Look at this script to see the build/run commands
    source run_tests.sh > outfile.log

As above, you can count the PASS/FAIL reports in outfile.log, and eyeball it.

    # Count the passes
    grep "PASS" outfile.log | wc -w
    # Count the failures (should be 0)
    grep "FAIL" outfile.log | wc -w

Currently the parallel tests print out other information (e.g. nesting
metadata) in outfile.log. That's OK -- but there should not be any FAILs or
obvious error messages.

# Step 3: Run the validation tests

Open a terminal in the [./tests/validation_tests](./tests/validation_tests)
folder, and run:

    Rscript run_validations.R

This will run over a dozen tests in examples/, and report one or more PASS
statements for each. There should be no FAIL statements if your install is
working. If you are missing some prerequisites (e.g. for coarrays, or R
packages), then there will be some failures. 

The code also generates various figures in the relevant directories in
examples/, and these should be visually inspected to better understand the test
results. 

Note the validation test PASS/FAIL statements give a relatively crude indicator of
model performance. They generally only check a few aspects of the results, and
the PASS/FAIL criteria depend on arbitrary thresholds. Such tests are useful to
catch changes in the model behaviour, and catastrophic errors. However they are
not a replacement for eyeballing the figures and thinking about the results. 

# Running models
----------------

# Examples to consider

The validation tests provide templates for developing other models. They illustrate driver programs, compilation, and model post-processing. Some of the more practical examples include:

* [./examples/nthmp/BP09/](./examples/nthmp/BP09) which simulates the Okushiri Island tsunami using multiple nested grids, and compares with observations
* [./examples/nthmp/Tauranga_harbour_Tohoku_tsunami/](./examples/nthmp/Tauranga_harbour_Tohoku_tsunami) which simulates the Tohoku tsunami at Tauranga harbour, NZ, and compares with velocity and tide-gauge observations. 
* [./examples/periodic_multidomain/](./examples/periodic_multidomain) which illustrates a global multidomain with periodic east-west boundaries

The above models can be run with OpenMP and/or coarrays, and illustrate use of the multidomain class. Another useful example is:

* [./examples/generic_model/](./examples/generic_model) which can run basic single-grid spherical coordinate models, e.g. for oceanic-scale tsunami modelling. This supports OpenMP but not coarrays.

# Advanced compilation

The script [./src/src_standard_compiler_var](./src/src_standard_compiler_var) is included in build scripts to define compiler variables (e.g. see the example models above). The user can override most variables in this file by pre-defining them before it is included. This is required, e.g. to use compilers other than gfortran.

A number of preprocessor variables can be defined to control features of the code. These are controlled by defining the variable `SWALS_PREPROCESSOR_FLAGS` in the makefile. See makefiles in the example projects for illustrations of their use. The most important cases include

- `-DTIMER` Time sections of the code and report on how long they take, generally in the multidomain log file. This is useful for understanding run-times. It is also used for static load-balancing calculations - see the function `make_load_balance_partition` in [./plot.R](./plot.R) which uses multidomain log-files to compute a more balanced distribution of work for multi-image coarray runs.
- `-DSPHERICAL` Assume spherical coordinates. Otherwise cartesian coordinates are used
- `-DREALFLOAT` Use single precision for all reals. Otherwise double-precision is used. Beware the nonlinear solvers generally need double-precision for accuracy.
- `-DNOCORIOLIS` Do not include Coriolis terms in spherical coordinates. By default, Coriolis terms are used when `-DSPHERICAL` is defined. They can ONLY be used in conjunction with spherical coordinates. But even in this case, sometimes it is useful to turn them off, hence this variable.
- `-DCOARRAY` Build with coarray support (alternatively the coarray calls may be replaced with MPI using flags below). Without this flag, only OpenMP parallisation is used.
- `-DCOARRAY_USE_MPI_FOR_INTENSIVE_COMMS` This uses MPI instead of coarrays for communication, which is useful because compiler support for coarrays is quite variable. If using this option you MUST also use -DCOARRAY (yes, even though it replaces coarrays).
- `-DCOARRAY_PROVIDE_CO_ROUTINES` Provide implementations of coarray collectives using MPI. This is required for compilers that do not support Fortran coarray collectives such as co_min, co_max, co_sum, etc. If using this option you MUST also use -DCOARRAY.
- `-DLOCAL_TIMESTEPPING_PARTITIONED_DOMAINS` Allow nonlinear domains inside a multidomain to take larger timesteps than suggested by `domain%timestepping_refinement_factor`, if this would be stable according to their own cfl-limit. This can speed up model runs, but also introduces load imbalance. The load imbalance can be dealt with by providing a load_balance_partition file (e.g. `md%load_balance_file="load_balance_partition.txt"`), which can be generated from a preliminary model run. See `make_load_balance_partition` in [./plot.R](./plot.R).
- `-DNETCDF4_GRIDS` Use the HDF5-based netcdf4 format for grid-file output. This requires that the netcdf library is compiled with netcdf4 support -- if not it will cause compilation to fail.

Other options that are less often useful include:

- `-DNOFRICTION` Do not use friction terms in the nonlinear shallow water equations. This can improve the speed for frictionless cases.
- `-DNOOPENMP` Do not use the openmp library, not even for timing the code. In this case, the timer will report the CPU time for all cores, not the wallclock time. This can occasionally be useful if you must avoid using openmp.
- `-DNONETCDF` Do not use netcdf for model outputs. This is useful if you cannot build with netcdf for some reason. However the output format is poorly supported, this is really just for testing if you'd like to bypass netcdf troubles.
- `-DMPI` Use `MPI_Abort` in the generic stop routine. This can probably be removed?
- `-DDEBUG_ARRAY` Add an array `domain%debug_array(nx, ny)` to every domain, which is written to netcdf at every output timestep. This can provide a scratch space for debugging.
