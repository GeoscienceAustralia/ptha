SWALS
-----

SWALS is a light-weight linear shallow water equations solver. It can run in
Cartesian or spherical coordinates, and output stage and depth-integrated
velocities over time, either as grids, and/or at tidal gauges. It can also
output the peak stage for an entire run. 

SWALS can run either on a single core, or in parallel using openmp. In the
latter case, speedups of around 4x are common (say on 6 cores), but you should
not expect good performance gains by extending this to many CPUs. The authors'
use cases have mainly involved running many model scenarios, where good
performance in the single-core scenario is most relevant. However, when just a
few runs are of interest, the parallel speedup is convenient. 

The code can run using either single or double precision arithmetic (see
section below on preprocessing flags). Single precision is usually fine (as
well as being faster, and using less memory), but users employing single
precision should confirm it is sufficient for their application, by doing a few
comparisons with the double precision version.

Model run-times will of course depend on the hardware you use. For example,
SWALS has been used to run a Pacific wide tsunami model at 4-arc-minute
resolution for 24 hours, in around 220s (with no IO) on a single core of my 3
year old desktop (or on raijin on the NCI). Times are a bit longer with IO,
depending on how much output is generated. As another example, on a single core
of the same hardware, we regularly run 'nearly global' models (360 degrees
longitude, latitude ranging from -72 to +65) at 1-arc-minute resolution, with
run times being about 75% of the real time (i.e.  around 18 hours to simulate
24 hours of tsunami propagation). We often simultaneously run 16 models like the
latter on a single node of the NCI (16 cores, 64GB total RAM) with similar run
times. In this case, note that using the 'numactl' program to enforce
reasonable CPU sharing can lead to better load-balancing and faster run times.


Testing
-------

SWALS includes a unit test suite in [test/unit_tests](test/unit_tests). These
are most useful for ensuring that your install is working, and to help us avoid
accidentally breaking features as the code evolves. When trying to install SWALS,
the first thing to do is try to compile and successfully run the unit tests.

It also ships with some analytical tests, including: 
* Three 1D solutions of wave shoaling over a piecewise linear beach, from the
 'national tsunami hazard mapping program' test suite
[(examples/nthmp/BP02)](examples/nthmp/BP02). To my knowledge, these are the
only ones from that test suite which are suitable for purely linear shallow
water equations solvers. The solution was obtained from:
https://github.com/rjleveque/nthmp-benchmark-problems .
* A 2D analytical solution of wave scattering around a conical island
[(examples/circular_island)](examples/circular_island). Unlike
the above tests, this test is 2D situation with variable topography. 

Mass conservation tracking is also implemented (both the volume in the domain,
and the time-integrated fluxes through the boundaries). The example scripts
show how to print out the difference between these quantities at regular
intervals (any change from a constant implies mass conservation errors). You
should not expect any mass conservation errors, (although tiny, non-systematic
floating-point related changes in the mass balance are normal).

A range of other tests are currently not provided here to avoid the need to
distribute large elevation data. It has been compared with the linear solver in
JAGURS on a realistic spherical coordinates case, giving identical answers up
to floating point precision (which it should, since the same widely used leap-frog 
algorithm is implemented in both). 

As with all PDE based numerical modelling, we strongly encourage you to do
convergence testing on all applications, to check the sensitivity of results to
grid size.


Getting started
---------------

SWALS requires existing installs of gfortran (>= 4.7), gdal, and netcdf, with
the latter compiled with the same version of gfortran that you use to compile
SWALS. It is setup to run on linux type environments. 

If you have those dependencies, then try to compile and run the unit tests in
[test/unit_tests](test/unit_tests). The programs therein should print many
'PASS' statements, and no 'FAIL' statements or other errors. 

If that works, then try running the test cases in
[examples/nthmp/BP02](examples/nthmp/BP02) and
[examples/circular_island](examples/circular_island), which will compare SWALS
with analytical solutions (and experimental data) for well known
linear-shallow-water test problem. Plotting requires that R is installed on
your system, along with the ncdf4 R package, although you can run the models
alone without R.

For spherical coordinate applications with an internally generated tsunami
source, see the example in [examples/generic_model](examples/circular_island).
Currently we do not provide the data to run that example (since it would be large), 
but it can be run if the user provide raster files with the elevation and initial stage.


Compilation details
-------------------
A number of preprocessor options can be provided to the compiler to control features of the code.
See makefiles in the example projects for illustrations of their use

    -DREALFLOAT (use single precision for all reals. Otherwise double-precision is used)
    -DSPHERICAL (assume spherical coordinates. Otherwise cartesian coordinates are used)
    -DNONETCDF (do not use netcdf for tide gauge outputs. This is useful if you cannot build with netcdf for some reason.)
    -DNOOPENMP (do not use the openmp library, not even for timing the code. In this case, the timer will report the CPU time for all cores, not the wallclock time. This can occasionally be useful if you must avoid using openmp.)
    -DTIMER (time sections of the code and report on how long they take)


If compiling on NCI, see the additional information in
[test/unit_tests](test/unit_tests), which is equally applicable to the other
examples.
