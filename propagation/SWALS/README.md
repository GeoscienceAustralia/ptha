SWALS
-----

SWALS includes a number of "Shallow WAter Like Solvers". These solve the linear
and nonlinear shallow water equations (the latter with Manning friction) using
either Cartesian or Spherical Coordinates. 

The shallow water equations are widely used to model large scale hydrodynamic
processes (such as tsunamis, floods), for which the horizontal length-scales of
interest are large compared with the flow depth. 

Features
--------

# User provides the driver program
SWALS requires the user to provide a high-level Fortran driver program to
specify the model domain(s), initial and boundary conditions, and IO. 

This makes model setup quite flexible.  

# Test suites 
SWALS includes a validation test suite with a variety of (mostly tsunami-type)
problems. In addition to illustrating the kinds of problems that can be solved,
the validation test suite serves to illustrate approaches to writing the driver
programs, compiling the code, and working with outputs. 

The code also has a unit-test suite, and a separate parallel unit-test suite.
These test suites are useful to confirm that the installation is working - and
reduce the risk of updates breaking the code.

# Two-way nesting
Two-way nesting is supported using a "multidomain" class which contains one or
more single-grid "domains". 

Where domains with different resolution overlap, the finer one is assumed to be
the "priority domain" which represents the real flow. Domains with the same
resolution should not overlap (because the code would not know which one should
be the priority domain), but their boundaries can touch. 

The code automatically generates the domain halos and data structures required
to populate halos using information from the corresponding priority domain. The
nested domains can use different solvers (e.g. linear or various non-linear
algorithms), and they can take different time-steps.

No adaptive-mesh refinement is used (see GEOCLAW and Basilisk for well-known
models which do that). However, if some domains are known to be inactive for
part of the simulation, they can be "switched off" (i.e. no evolution of the
flow) before a specified time. 

Flux-correction is used to ensure mass conservation between nested grids, and
mass tracking routines are provided to confirm that the total mass in the
multidomain is consistent with fluxes through the multidomain boundaries.

# Multiple types of solvers
The domains contained within a multidomain can use different solvers, and take
different timesteps. In general the domain solver algorithm should be chosen
depending on the problem.

The linear solver is well suited to modelling oceanic-scale
earthquake-generated tsunamis, in regions where the wave amplitude is small
compared with the depth. However the linear equations are inappropriate for
inundation or nearshore tsunami modelling, where non-linearity matters.

Conversely, the nonlinear solvers mostly employ 2nd-order shock-capturing
finite-volume methods, along the lines of those implemented in the
unstructured-mesh ANUGA code (although SWALS uses structured grids). These
methods are quite robust, but inefficient for modelling deep ocean tsunamis at
global-type spatial scales when a linear solver would suffice (and would be
much faster -- say by "several factors of 10"). 

# Parallel support with OpenMP (shared memory) and/or Fortran coarrays (distributed or shared memory)
SWALS can be run either on a single core, or in parallel using OpenMP and/or
Fortran coarrays. 

We have observed good strong-scaling both with on-node (6-16 cores) and
distributed memory parallelism over a nodes (up to 16 nodes = 256 cores on
Raijin). There is nothing preventing larger runs. However the capacity of any
model to efficiently use many cores is heavily dependent on the problem size,
load balancing, and the relative use of OpenMP threads and coarrays -- factors
within the users control. Thus the model scaling will vary in a case-by-case
manner.

The use of more coarray images leads to relatively more 'halos cells' in the
domain, and with enough coarray images this becomes inefficient. The OpenMP
parallelization does not require halos, but involves more on-node
synchronization which is intrinsically less efficient. The optimal balance
between OpenMP threads and coarray images involves balancing these two sources
of inefficiency. In the authors experience, it is often better to just use 2-4
OpenMP threads per coarray image (rather than 1 coarray image per node). But
this will be problem and hardware dependent. 

The test-suite includes a few comparisons of models using OpenMP and mixed
OpenMP-coarrays. Because floating point operations are not associative, some
numerical differences can arise using different parallel approaches, and
ultimately this may lead to differences in the results depending on the
parallel approach. However we do not expect substantial differences - results
should be "the same" or "similar in all practical respects". Quantities most
likely to show differences are those that are sensitive to round-off type
changes in variables. For example mass conservation errors will ideally be
negligible (i.e. completely reflecting floating point limitations) - so these
might vary among parallel runs, which is fine so long as they remain negligibly
small. As another example, SWALS treats a cell as wet or dry depending on
whether its flow depth exceeds a threshold (1.0e-05 m). This wet/dry status
could be sensitive to round-off in extremely shallow cells, and sometimes this
shows up as differences in reported runup maxima over time.


# Double or single precision
SWALS can run using either double precision (recommended) or single precision
reals. For many problems single precision is not sufficient - but this can be
checked on a case-by-case basis by comparing with the double-precision results.
In general, problems using nesting and/or the nonlinear solvers over large
depths tend to require double-precision.

# Gauge and Grid output
SWALS can output stage and depth-integrated velocities over time, as grids
and/or gauges, in netcdf format. It can also output the peak stage for an
entire run.

Installation prerequisites
--------------------------

## Fortran compiler
Compilation requires a version of GNU-Fortran, combined with OpenCoarrays to
provide distributed memory coarrary support. The author hasn't comprehensively
tested various compiler versions, but it is likely that gfortran versions 5 or
greater will work (without coarrays), but OpenCoarrays may require a more
recent gfortran.

The code has also been compiled and run with Intel-Fortran 19, which has some
support for coarrays. However this requires different preprocessor options to
run efficiently, see below. It also requires adjustments to the build script.

## C compiler
The code also makes limited use of the C-language to interface with the GDAL
library (thus providing flexible raster-based input), so a C compiler is
required (typically gcc, or use Intel's icc with Intel-Fortran).

## Make
We also require the "make" program for compilation.

## NetCDF and GDAL
"NetCDF" and "GDAL" must also be installed. The command-line programs
"gdal-config" and "nc-config" or "nf-config" should be available, because the
makefiles use them to identify the associated libraries. The netcdf install
should include fortran support, and this requires building netcdf with the same
fortran compiler as you use to build SWALS.

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
    install.packages('netcdf4', 'fields', 'raster', 'sp', 'rgdal', 'maptools', 'rasterVis')


Compiling and testing
---------------------

The standard build scripts are based on src/src_standard_compiler_var and
src/src_make_commands. These scripts are 'included' in application specific
makefiles.  

Variables in src_standard_compiler_var can be overridden by defining them
before the include. This is required e.g. to use a different compiler or
non-standard library locations. 

See the validation test suite for examples.

# Step 1: Run the unit-tests

The unit test suite is in [test/unit_tests](test/unit_tests). These are most
useful for ensuring that your install is working, and to help us avoid
accidentally breaking features as the code evolves. When trying to install
SWALS, the first thing to do is try to compile and successfully run the unit
tests. To do this, open a terminal in that directory and do

    make -B -f make_test
    ./unit_tests > outfile.log

This should take a minute, and write a lot of "PASS" statements to outfile.log. There
should not be any "FAIL" statements or other warning messages. Open up outfile.log to check.
You can count PASS or FAIL statments like this:

    # Count the passes
    grep "PASS" outfile.log | wc -w
    # Count the failures (should be 0)
    grep "FAIL" outfile.log | wc -w

# Step 2: Run the parallel unit-tests

The parallel unit-test suite is in [test/parallel_tests](test/parallel_tests). Note
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

Open a terminal in the tests/validation_tests folder, and run:

    Rscript run_validations.R

This will run over a dozen tests in examples/, and report one or more PASS
statements for each. There should be no FAIL statements if your install is
working. If you are missing some prerequisites (e.g. for coarrays), then there
will be some failures. 

The code also generates various figures in the relevant directories in
examples/, and these should be visually inspected to better understand the test
results. 

The validation test PASS/FAIL statements give a relatively crude indicator of
model performance. They generally only check a few aspects of the results, and
the PASS/FAIL criteria depend on arbitrary thresholds. Such tests are useful to
catch changes in the model behaviour, and catastrophic errors. However they are
not a replacement for eyeballing the figures and thinking about the results. 

# Running models
----------------

# Examples to consider

The validation tests provide useful templates for developing other models. They illustrate driver programs, compilation, and model post-processing. Particularly useful examples include:

* ./examples/nthmp/BP09/ which simulates the Okushiri Island tsunami using multiple nested grids, and compares with observations
* ./examples/nthmp/Tauranga_harbour_Tohoku_tsunami/ which simulates the Tohoku tsunami at Tauranga harbour, NZ, and compares with velocity and tide-gauge observations. 
* ./examples/periodic_multidomain/ which illustrates a global multidomain with periodic east-west boundaries

The above models can be run with OpenMP and/or coarrays, and illustrate use of the multidomain class. Another useful example is:

* ./examples/generic_model/ which can run basic single-grid spherical coordinate models, e.g. for oceanic-scale tsunami modelling. This supports OpenMP but not coarrays.

# Advanced compilation

The script src/src_standard_compiler_var is included in build scripts to define compiler variables (e.g. see the example models above). The user can override most variables in this file by pre-defining them before it is included. This is required, e.g. to use compilers other than gfortran.

A number of preprocessor variables can be defined to control features of the code. These are controlled by defining the variable SWALS_PREPROCESSOR_VAR in the makefile. See makefiles in the example projects for illustrations of their use. The most important cases include

    -DSPHERICAL (assume spherical coordinates. Otherwise cartesian coordinates are used)
    -DTIMER (time sections of the code and report on how long they take)
    -DCOARRAY (build with coarray support. Otherwise only OpenMP parallisation can be used)
    -DREALFLOAT (use single precision for all reals. Otherwise double-precision is used)
    -DNOCORIOLIS (do not include Coriolis terms in spherical coordinates. By default, Coriolis terms are used when -DSPHERICAL is defined. They can ONLY be used in conjunction with spherical coordinates. But even in this case, sometimes it is useful to turn them off, hence this variable)
    -DCOARRAY_USE_MPI_FOR_INTENSIVE_COMMS (this uses MPI instead of coarrays for communication in the main evolve loop. It improved performance using Intel-Fortran 2019, which does not have great coarray performance)
    -DCOARRAY_PROVIDE_CO_ROUTINES (provide implementations of coarray collectives using MPI. This is required for Intel-Fortran 2019, which does not support Fortran coarray collectives such as co_min, co_max, co_sum, etc)

Other options that are less often useful include:

    -DNOFRICTION (do not use friction terms in the nonlinear shallow water equations. This can improve the speed for frictionless cases)
    -DNOOPENMP (do not use the openmp library, not even for timing the code. In this case, the timer will report the CPU time for all cores, not the wallclock time. This can occasionally be useful if you must avoid using openmp.)
    -DNONETCDF (do not use netcdf for tide gauge outputs. This is useful if you cannot build with netcdf for some reason.)
    -DMPI (use MPI_Abort in the generic stop routine. This can probably be removed?)

