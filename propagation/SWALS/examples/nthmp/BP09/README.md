# NTHMP Benchmark problem 9: Okushiri tsunami simulation and comparison with field data

This problem simulates the 12 July 1993 Hokkaido-Nansei-Oki tsunami (often referred to as the Okushiri tsunami), and compare the model with field observations of tsunami runup.

This test problem is from the NTHMP benchmark suite. The test data and a problem description is available in [Randy LeVeque's repository](https://github.com/rjleveque/nthmp-benchmark-problems/tree/master/BP09-FrankG-Okushiri_island/). 

The [SWALS_model](BP09.f90) simulates this problem for one hour following the earthquake. The multidomain uses 6 nested domains and spherical coordinates. The outer domain uses the `linear` solver, while all nested domains use the default nonlinear solver `rk2`. 

## Issues with the benchmark problem datasets

Although models reasonably produce observations for this benchmark problem, the underlying data also has some weaknesses (also discussed in the [GEOCLAW benchmark report](https://depts.washington.edu/clawpack/links/nthmp-benchmarks/geoclaw-results.pdf):
* The elevation data has clear discontinuities and mismatches between neighbouring grids. For the SWALS test, some effort was made to preprocess DEMs to reduce discontinuities, but elevation artefacts definitely remain and must be affecting the model.
* The runup observations have substantial location errors, relative to the DEMs. To partially correct for this we have estimated offsets for each dataset, which are applied in [plot.R](plot.R). This improves the positions, but does not resolve all the location errors.

While this remains a very useful benchmark problem, caution should be applied to avoid over-interpreting results. 

## Stability of the results with different domain partitions and timestepping

As part of this test problem the [run_model.sh](run_model.sh) script re-runs the problem with several different setups.
1. The first run uses openmp alone
2. The second run uses a mixture of openmp and MPI, with default partitioning of domains.
3. The third run also uses a mixture of openmp and MPI, with prescribed partitioning of domains, and allows nonlinear domains to take longer timesteps if their CFL allows it. This option is enabled by compiling SWALS with `-DLOCAL_TIMESTEP_PARTIONED_DOMAINS`.

The test code checks that all of these models give similar results. While we can force the results of runs 1 and 2 to be identical (discussed below) with the default setup we expect small differences, because:
* For runs 1 and 2 we do not specify how the domains should be partitioned. SWALS default behaviour is for run 1 to use the original 6 domains, while run 2 splits each domain into as many pieces as there are MPI ranks and spreads them among processes. Due to the use of finite precision floating point numbers to define the domains, the partitioned domain coordinates are not __exactly__ equivalent to the original domain coordinates (although the differences are tiny). Thus, small differences in the solutions occur. 
    * It is possible to force runs 1 and 2 to to give __identical__ results by prescribing the same domain partition in both cases (i.e. setting `md%load_balance_file`). This is good practice in applications. 
      * Such a test is implemented in an alternative script [run_model_exact_reproduce.sh](run_model_exact_reproduce.sh). For brevity that is not run by the automated test suite, but an analogous test is part of the [paraboloid_bowl](../../paraboloid_bowl) tests.
* For run 3 we specify how the domain should be partitioned, and additionally allow the nonlinear domains to take longer time-steps if stability permits. (This is not done for leapfrog timestepping methods, because theoretically they rely on a fixed timestep). As a result there are also small differences with the previous runs. Again the test code confirms that the differences are very small.

Note that while we can force runs 1 and 2 to give the same results by prescribing domain partition, is not possible for force run 3 to give identical results as either run 1 or run 2. This is because by allowing different timestepping in run 3, the numerical method has changed. 

Another interesting test (not run by default) is to repeat run 1 while adding a tiny perturbation to the model elevations (1e-10 m). This leads to differences comparable to those between runs 1 and 2, and further highlights that tiny floating point differences in the initial conditions lead to small differences in the results.

