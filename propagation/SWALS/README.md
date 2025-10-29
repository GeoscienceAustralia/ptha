# SWALS
-----

Shallow WAter Like Solvers (SWALS) computes solutions to several [variants of
the 2D shallow-water equations](./SOLVERS.md) (linear/nonlinear/dispersive) in Cartesian
and Spherical coordinates, on "multidomains" represented as a set of connected and/or
nested rectangular grid domains. 

A number of [different numerical methods](./SOLVERS.md) are implemented,
suitable for a range of flow regimes, with particular emphasis on tsunami-like
problems. This includes: 
* [shock-capturing finite volume schemes](https://github.com/GeoscienceAustralia/ptha/blob/master/propagation/SWALS/SOLVERS.md#the-finite-volume-solvers)
* [classical leapfrog schemes](https://github.com/GeoscienceAustralia/ptha/blob/master/propagation/SWALS/SOLVERS.md#the-leapfrog-schemes ),
* The [CLIFFS scheme](https://github.com/GeoscienceAustralia/ptha/blob/master/propagation/SWALS/SOLVERS.md#the-cliffs-solver)
which was [developed by Elena Tolkova](https://github.com/Delta-function/cliffs-src). CLIFFS is similar to the well-known MOST tsunami solver, but uses a different wetting and drying scheme. 

The figure below shows an example SWALS model, which uses several levels of  nested grids to simulate tsunamis at tide gauges in southeast and west Australia (see [here](https://doi.org/10.1029/2025JB031949)). 

![SWALS model used to simulate tsunamis at tide gauges in southeast and west Australia.](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/Davies2025_Australia_wide_model_structure.png)

A SWALS model with a more complex nesting structure is shown below. This was used to simulate tsunami inundation throughout NSW in Australia ([see here](https://www.researchgate.net/publication/396141903_NSW-wide_probabilistic_tsunami_inundation_hazards_from_subduction_earthquakes)).

![SWALS model used to simulate tsunami inundation throughout NSW in Australia.](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/Structure-of-the-NSW-wide-inundation-model-A-The-global-extent-with-large-boxes_W640.jpg)

Two-way nesting allows for the use of higher-resolution domains in specified
areas, without hard limits on the number of domains or depth of refinement.
Nested domains may use different numerical solvers, take different timesteps,
and have grid sizes that are any integer divisor of the coarser domain(s) they
communicate with. 
* The finest-resolution domain at any particular point is the "priority domain" and contains the SWALS numerical solution. This is communicated to "halo regions" in neighbouring domains to enable seamless evolution of the flow.
* Flux correction is used to enforce the conservation of mass and advected momentum through nested domain interfaces (for schemes that would conserve these quantities on a single grid). In practice we obtain excellent mass conservation in complex multidomains when using a combination of leapfrog and finite-volume schemes. 

Parallel computation (shared and distributed memory CPU) is supported with a
mixture of MPI (or Fortran coarrays) and openmp. Domains can be automatically split
between MPI ranks, or the partition can be specified by the user at run-time.
Static load balancing can be used to improve the efficiency of large parallel
jobs.

The code includes various test suits that [can be run automatically](#compiling-and-testing), including a
[unit test suite](./tests/unit_tests), a [parallel unit test suite](./tests/parallel_tests),
 and a [validation test suite](./tests/validation_tests). The latter focusses 
on tsunami type problems; see [here for various NTHMP tests](./examples/nthmp) 
(which are well known in the tsunami community), [here](./examples) for other problems, and [here](./examples/dispersive/) for tests using dispersion. 

Publications using SWALS include:
* [A paper](https://www.frontiersin.org/articles/10.3389/feart.2020.598235/full) using SWALS to model 5 historic tsunamis in Australia.
It includes discussion of the range of equations that are solved, the energy conservation properties of different solvers, and their degree of agreement with observations. 
* [A paper](https://doi.org/10.1029/2025JB031949) using SWALS to model 14 historical tsunamis at 32 tide gauges in southeast and west Australia. Observed tsunamis are compared with models using both published source inversions and random PTHA18 scenarios.
* [A paper](https://doi.org/10.1093/gji/ggac140) using SWALS model inundation hazards in Tongatapu. It includes comparison with observations of 5 historical tsunamis at Nuku'alofa tide gauge (using PTHA18 scenarios as the source models).
* [A technical report](http://dx.doi.org/10.26186/150015) and [conference paper](https://icce-ojs-tamu.tdl.org/icce/article/view/12657/11930) using SWALS to model tsunamis hazards in Western Australia (Geraldton to Dunsborough). The technical report includes a comparison of models with observations of the 2004 and 2005 tsunamis at many tide gauges, and some model intercomparisons.
* The linear solver in SWALS was used extensively for the [2018 Australian Probabilistic Tsunami Hazard Assessment](http://dx.doi.org/10.11636/Record.2018.041) and two associated papers: [this one in GJI](https://doi.org/10.1093/gji/ggz260) and  [this one in PAGEOPH](https://link.springer.com/article/10.1007/s00024-019-02299-w). 
* [A conference paper](https://www.researchgate.net/publication/396141903_NSW-wide_probabilistic_tsunami_inundation_hazards_from_subduction_earthquakes) using SWALS to model tsunami hazards in NSW, including a comparison of the model with observations of the 2021 Kermadec tsunami.
* [A conference paper](https://www.researchgate.net/publication/396142067_Accommodating_Spatially_Varying_Tidal_Planes_in_Tsunami_Hazard_Assessments) using SWALS to model tsunami hazards near Gladstone in QLD, including a comparison of the model with observations of the 2007 Solomon tsunami.

## Installation prerequisites
--------------------------

## Fortran compiler
Compilation requires a version of GNU-Fortran or ifort, although the test-suite
is setup for gfortran by default. The author hasn't comprehensively tested
various compiler versions, but as of 2025 it is being used with versions of
`gfortran`, `ifort` and `ifx` released between 2021-2025. 

The code can run in parallel with both shared memory (openmp) and distributed
memory approaches. Distributed parallel support requires that MPI is installed,
and/or that the compiler has adequate coarray support. In general all these
parallel approaches can be used together.

## C compiler
The code also makes limited use of C to interface with the GDAL library (thus
providing flexible raster-based input). Thus a C compiler is required,
typically gcc.

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
    install.packages(c('ncdf4', 'fields', 'raster', 'sp'))


# Compiling and testing
---------------------

## Step 0: Adapting the build scripts for your machine 

The standard build scripts are based on
[src_standard_compiler_var](./src/src_standard_compiler_var) and
[src_make_commands](./src/src_make_commands). These are 'included' in
the application specific makefiles which compile SWALS.

The script [src_standard_compiler_var](./src/src_standard_compiler_var) tells SWALS
where to find your compilers and other libraries. To make it easier to move between machines,
[src_standard_compiler_var](./src/src_standard_compiler_var) simply
includes another machine specific script. A few variants of these are provided, which use
[gfortran](./src/src_standard_compiler_var_gfortran) and
[ifort](./src/src_standard_compiler_var_ifort), and one that works on 
[the gadi machine on NCI](./src/src_standard_compiler_var_NCI_gadi_ifort).
Make sure you edit [src_standard_compiler_var](./src/src_standard_compiler_var)
so the variables are suitable for your machine.

Variables in [src_standard_compiler_var](./src/src_standard_compiler_var) can
be overridden in particular applications, by defining them in the
application-specific makefile (see the examples). This is required in many
situations (e.g. to use spherical coordinates; different compiler
options; non-standard library locations). There are many examples of
application makefiles in the [examples folder](./examples/) (look for files
with names beginning with make\_). For instance [this
makefile](./examples/nthmp/BP09/make_BP09_coarray) builds a
spherical-coordinates model with distributed-memory-parallel support, while
[this makefile](examples/circular_island/make_circular_island) builds an
openmp-only model with Cartesian coordinates using single-precision reals and
code-timing. To make a new makefile, you often just need to copy an existing
one and change the 'mymodel' variable to correspond to your main f90 file.

See the validation test suite for examples, and documentation of compiler options below.

## Step 1: Run the unit-tests

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

## Step 2: Run the parallel unit-tests

The parallel unit-test suite is in [./tests/parallel_tests](./tests/parallel_tests). Note
this only tests the distributed memory parallel code (MPI or coarray), not
OpenMP; both are tested in the validation suite. 

To run the tests, open a terminal in that directory and do:

    # Look at this script to see the build/run commands
    source run_test.sh > outfile.log

As above, you can count the PASS/FAIL reports in outfile.log, and eyeball it.

    # Count the passes
    grep "PASS" outfile.log | wc -w
    # Count the failures (should be 0)
    grep "FAIL" outfile.log | wc -w

Currently the parallel tests print out other information (e.g. nesting
metadata) in outfile.log. That's OK -- but there should not be any FAILs or
obvious error messages.

## Step 3: Run the validation tests
These tests include a range of analytical flow solutions, and comparisons with laboratory and field data. 
They also check that results are consistent when running in parallel in different ways.

Before running the validation tests, you might need to modify
[./src/test_run_commands](./src/test_run_commands) to tell the validation test suite how to run MPI/openmp
on your machine. This file points to another script to make it easier to move between machines.
For instance the script [./src/test_run_commands_basic](./src/test_run_commands_basic) works on
my home machine, while [./src/test_run_commands_gadi_1node](./src/test_run_commands_gadi_1node) works
on a single node of NCI's gadi machine.

To run the validation tests, open a terminal in the
[./tests/validation_tests](./tests/validation_tests) folder, and run:

    Rscript run_validations.R

This will run over a dozen tests in examples/, and report one or more PASS
statements for each. There should be no FAIL statements if your install is
working. If you are missing some prerequisites (e.g. R packages), then there
will be some failures. 

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

## Examples to consider

A simple example that illustrates basic useage with a multidomain is in [./examples/radial_dam_break](./examples/radial_dam_break).

The validation tests  in [./examples/](./examples) provide templates for developing other models. They illustrate driver programs, compilation, and model post-processing. Some of the more practical examples include:

* [./examples/nthmp/BP09/](./examples/nthmp/BP09) which simulates the Okushiri Island tsunami using multiple nested grids, and compares with observations
* [./examples/nthmp/Tauranga_harbour_Tohoku_tsunami/](./examples/nthmp/Tauranga_harbour_Tohoku_tsunami) which simulates the Tohoku tsunami at Tauranga harbour, NZ, and compares with velocity and tide-gauge observations. 
* [./examples/periodic_multidomain/](./examples/periodic_multidomain) which illustrates a global multidomain with periodic east-west boundaries. This also optionally permits a rise-time to be used in the earthquake co-seismic deformation.
* The model code used in [this study modelling historic tsunamis in Australia](https://www.frontiersin.org/articles/10.3389/feart.2020.598235/full) can be [found here](../../misc/nearshore_testing_2020/). It implements a global-to-local scale model with an initial earthquake forcing.
* A model used for [work in Western Australia](https://icce-ojs-tamu.tdl.org/icce/article/view/12657/11930) that makes complicated nesting relatively easy can be [found here](https://github.com/GeoscienceAustralia/ptha/tree/master/misc/SW_WA_2021_2024/greater_perth_revised2023).

The above models can be run with OpenMP and/or MPI (or coarrays), and illustrate use of the multidomain class. 

Another useful example is:

* [./examples/generic_model/](./examples/generic_model) which can run basic single-grid spherical coordinate models, e.g. for oceanic-scale tsunami modelling. It permits a rise-time to be applied to the initial co-seismic deformation. This supports OpenMP but not MPI/coarray (for a broadly similar model with distributed memory support, see [./examples/periodic_multidomain/](./examples/periodic_multidomain)). 


## Advanced compilation

The script [./src/src_standard_compiler_var](./src/src_standard_compiler_var) is included in build scripts to define compiler variables (e.g. see the example models above). The user can override most variables in this file by pre-defining them before it is included. This is required, e.g. to use spherical coordinates, single precision, unusual compiler options, etc.

A number of preprocessor variables can be defined to control features of the code. These are controlled by defining the variable `SWALS_PREPROCESSOR_FLAGS` in the makefile. See makefiles in the example projects for illustrations of their use. The most important cases include

- `-DTIMER` Time sections of the code and report on how long they take, generally in the multidomain log file. This is useful for understanding run-times. It is also used for static load-balancing calculations - see the function `make_load_balance_partition` in [./plot.R](./plot.R) which uses multidomain log-files to compute a more balanced distribution of work for multi-image MPI/coarray runs.
- `-DSPHERICAL` Assume spherical coordinates. Otherwise cartesian coordinates are used
- `-DREALFLOAT` Use single precision for all reals. Otherwise double-precision is used. This generally makes the code run faster, but beware - the nonlinear solvers generally need double-precision for accuracy, as do nested-grid models. If in doubt, test!
- `-DNOCORIOLIS` Do not include Coriolis terms in spherical coordinates. By default, Coriolis terms are used when `-DSPHERICAL` is defined (except in the CLIFFS solver). They can ONLY be used in conjunction with spherical coordinates. But even in this case, sometimes it is useful to turn them off, hence this variable.
- `-DCOARRAY` Build with distributed memory parallel support (in addition to shared memory support with openmp - which is enabled by default in any case). If ONLY this flag is provided then coarrays are used - alternatively, the coarray calls may be replaced with MPI by ALSO using flags below.
- `-DCOARRAY_USE_MPI_FOR_INTENSIVE_COMMS` This uses MPI instead of coarrays for communication, which is useful because compiler support for coarrays is quite variable, while MPI is very mature. If using this option you MUST also use `-DCOARRAY` (even though MPI replaces coarrays). 
- `-DCOARRAY_PROVIDE_CO_ROUTINES` Provide implementations of coarray collectives using MPI. This is required for compilers that do not support Fortran coarray collectives such as `co_min`, `co_max`, `co_sum`, etc. If using this option you MUST also use `-DCOARRAY`, and `-DCOARRAY_USE_MPI_FOR_INTENSIVE_COMMS` if using MPI. To be clear; **if this option and the previous 2 are provided**, then the compiler does NOT need to support coarrays (i.e. **all distributed-memory communication is done with MPI**).
- `-DLOCAL_TIMESTEPPING_PARTITIONED_DOMAINS` Allow domains in a multidomain to take larger timesteps than suggested by `domain%timestepping_refinement_factor`, if this would be stable according to their own cfl-limit. Currently it is only allowed for non-staggered-grid schemes (finite-volume and CLIFFS) because the logic of the staggered-grid schemes requires a fixed timestep. The use of local timestepping can speed up model runs, but also introduces load imbalance. The load imbalance can be dealt with by providing a load_balance_partition file (e.g. `md%load_balance_file="load_balance_partition.txt"`), which can be generated from a preliminary model run. See `make_load_balance_partition` in [./plot.R](./plot.R).
- `-DNETCDF4_GRIDS` Use the HDF5-based netcdf4 format for grid-file output. This requires that the netcdf library is compiled with netcdf4 support -- if not it will cause compilation to fail.
- `-DTRACK_MULTIDOMAIN_STABILITY` Insert calls to `check_multidomain_stability` into the multidomain evolve loop. This checks every domain in the multidomain for `NaNs` or unusually small/large values, repeatedly during each time-step. If issues are detected then the multidomain writes various output and calls `error-stop`. In practice this is useful to identify the root-cause of stability problems in complex models. Note that if instability is detected on one multidomain in a distributed parallel model, then the call to `error-stop` can cause stability issues to be detected on other images right after the call to the communication routines (tagged with `step-after-comms`). In this case, the multidomain that "really" caused the issue can be identified because its instability is not triggered right after the communication step; it is more likely to have a tag like `inner` or `innerB`. See `evolve_multidomain_one_step` for the tags that are used.
- `-DDEBUG_ARRAY` Add an array `domain%debug_array(nx, ny)` to every domain, which is written to netcdf at every output timestep. This can provide a scratch space for debugging.
- `-DDISPERSIVE_PEREGRINE_IN_FLUX_FORM` If any domain uses dispersion (`domain%use_dispersion = .true.`) then this will cause the code to use a two-term dispersive model (like Peregrine's dispersion but written in terms of the flux rather than the speed). Otherwise we use the simplest 1-term dispersive model written in terms of the flux (like JAGURS and many other tsunami codes), which is the faster option. These two options give similar results in many test problems. When they differ the simpler model is sometimes better.

Other options that are less often useful include:

- `-DNOFRICTION` Do not use friction terms in the nonlinear shallow water equations. This can improve the speed for frictionless cases.
- `-DNOOPENMP` Do not use the openmp library, not even for timing the code. In this case the code will run on a single core (or a single core per-MPI process), and an alternative timer is used which measures time somewhat differently. This can occasionally be useful if you must avoid using an openmp library (e.g. in some debugging contexts). Note this is NOT required to run the code single-threaded; in that case just set `OMP_NUM_THREADS=1`.
- `-DNONETCDF` Do not use netcdf for model outputs. This is occasionally useful if you cannot build with netcdf for some reason and want to test something else; however the output format is poorly supported and the validation tests will not pass. While occasionally useful this is NOT a sustainable way to use the code.
- `-DCOARRAY_MPI_USE_ALLTOALLV` This has an effect if used in conjunction with `-DCOARRAY_USE_MPI_FOR_INTENSIVE_COMMS`. It makes the code use Mpi_AlltoAllv to do communication in the main timestepping routine, instead of point-to-point Mpi_isend/Mpi_ireceive calls. Typically this is slower, but I have seen some problems where it was faster.
- `-DLEGACY_CORIOLIS_PARAMETER` Previously SWALS use a value for the earth's angular frequency assuming one rotation per day. Although this matches several other tsunami codes, it is slightly incorrect so was updated (the previous value was about 0.3% too small, because the earth rotates once in a sidereal day, which is slightly less than 24 hours). This flag makes SWALS use the slightly-too-small value. In principle this might help to compare with a previous run, or another model that uses the same treatment, although in practice I haven't seen a case where it matters. 
- `-DEVOLVE_TIMER` This adds an addition timer to every domain object, which times different parts of code within its time-stepping routine. This level of timing granularity isn't needed in general, but may be useful in understanding the performance of the code. The timer results are written to a domain-specific file within the domain output directory, named `Evolve_timer_details.txt`.
- `-DOLD_PROCESS_DATA_TO_SEND_B4FEB22` Use an older version of nesting, implemented prior to Feb 2022 (changes `two_way_nesting_comms_mod::process_data_to_send`). The newer version seems better for communication when a fine staggered-grid solver receives from a coarser solver, and similar in other cases.
- `-DALWAYS_FLUX_CORRECT_DRY_CELLS` Try to obtain bitwise reproducibility of old flux correction behaviour. Only use this with purely non-dispersive models (it will make some dispersive models unstable). After dispersion was added to the code (Oct 2025), the flux correction calculations became more involved, albeit for non-dispersive domains the older approach is logically equivalent. Despite this, some compilers give floating-point level differences using the older and newer flux correction code, even in the non-dispersive case (e.g. `gfortran` with [BP09](examples/nthmp/BP09)). This is presumably because the older code is easier to optimise. This preprocessor variable applies the old flux correction logic, and can eliminate those floating-point level differences.

## Source-code html documentation
-------------------------

HTML documentation of the source-code can be generated using [ford](https://github.com/Fortran-FOSS-Programmers/ford). Assuming ford is installed you can create the documentation by running:

    ford documentation_ford.md 

in the current directory. This will make a folder `doc`. You can then browse the documentation by opening `./doc/index.html` in your web-browser. It may be useful to start by looking at the derived type `domain_type` which solves the shallow water equations on a single grid, as well as the `multidomain_type` which holds multiple domains and manages communication between them.


