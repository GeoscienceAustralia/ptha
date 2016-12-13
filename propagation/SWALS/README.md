SWALS
-----

SWALS is a light-weight linear shallow water equations solver. It can run in
cartesian or spherical coordinates, and output stage and depth-integrated
velocities over time, either as grids, and/or at tidal gauges. It can also
output the peak stage for an entire run. 

SWALS can run either on a single core, or in parallel using openmp. In the
latter case, in the authors experience speedups of around 4x are common (say on
6 cores), but you should not expect good performance gains by extending this to
many CPUs. 

Either single or double precision arithmetic is supported, and while single
precision is usually fine (as well as being faster, and using less memory),
users employing single precision should confirm it is sufficient for their
application by a few comparisons with the double precision version.

Run-times will of course depend on the hardware you use. For example, SWALS has
been used to run a Pacific wide tsunami model at 4-arc-minute resolution for 24
hours, in around 220-300s (depending on IO choices) on a single core of my 3
year old desktop (or on raijin on the NCI). As another example, on a single
core of the same hardware, we regularly run 'nearly global' models (360 degrees
longitude, latitude ranging from -72 to +65) at 1-arc-minute resolution, with
run times being abour 75% of the real time (i.e. around 18 hours to simulate 24
hours of tsunami propagation). We often simulaneously run 16 models like the
latter on a single node of the NCI (16 cores, 64GB total RAM) with similar run
times. In this case, note that using the 'numactl' program to enforce
reasonable CPU sharing can lead to better load-balancing and faster run times.


Testing
-------

SWALS includes a unit test suite in tests/unit_tests. It also ships with three
analytical tests (to my knowledge, these are the only ones from the national
tsunami hazard mapping program test suite which are suitable for purely linear
shallow water equations solvers, see
https://github.com/rjleveque/nthmp-benchmark-problems ).

Mass conservation tracking is also implemented (both the volume in the domain,
and the time-integrated fluxes through the boundaries). The example scripts
show how to print out the difference between these quantities at regular
intervals (any change from a constant implies mass conservation errors). You
should not expect any mass conservation errors, (although tiny, non-systematic
floating-point related changes in the mass balance are normal).

A range of other tests are currently not provided here to avoid the need to
distribute large elevation data. It has been compared with the linear solver in
JAGURS on a realistic spherical coordinates case, giving identical answers up
to floating point precision (which it should, since the same leap-frog algorithm
is implemented in both).


Getting started
---------------

SWALS requires existing installs of gfortran, gdal, and netcdf, with the latter
compiled with the same version of gfortran that you use to compile SWALS. It is
setup to run on linux type environments. 

If you have those, then try to compile and run the unit tests in tests/unit_tests.
The programs therein should print many 'PASS' statements, and no 'FAIL' statements
or other errors. 

If that works, then try running the test case in examples/nthmp/BP02, which will
compare SWALS with analytical solutions (and experimental data) for a well known
linear-shallow-water test problem. Plotting requires that R is installed on your 
system.

For spherical coordinate applications, see the example in examples/generic_model.
