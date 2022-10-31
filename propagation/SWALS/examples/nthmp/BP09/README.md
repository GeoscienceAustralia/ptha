# NTHMP Benchmark problem 9: Okushiri tsunami simulation and comparison with field data

This problem simulates the 12 July 1993 Hokkaido-Nansei-Oki tsunami (often referred to as the Okushiri tsunami), and compares the model with tsunami runup observations.

This test problem is from the NTHMP benchmark suite. The test data and problem description are available in [Randy LeVeque's repository](https://github.com/rjleveque/nthmp-benchmark-problems/tree/master/BP09-FrankG-Okushiri_island/). 

The [SWALS_model](BP09.f90) simulates the tsunami for one hour following the earthquake. The multidomain contains six domains in spherical coordinates. The outer domain uses the `linear` solver, while all other domains use the default nonlinear solver `rk2`. 

![Figure 1: Elevation data and multidomain design](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/nthmp/BP09/elevation_okushiri_lowresolution_omp.png)

## Issues with the benchmark problem datasets

For this benchmark we expect models to reproduce the observations quite well, but the underlying data has some weaknesses (also discussed in the [GEOCLAW benchmark report](https://depts.washington.edu/clawpack/links/nthmp-benchmarks/geoclaw-results.pdf)):

* The elevation has clear discontinuities and mismatches between neighbouring grids. For the SWALS test, some effort was made to preprocess DEMs to reduce mismatches, but elevation artefacts definitely remain.
* The runup observations have substantial location errors, relative to the DEMs. To partially correct for this we have estimated location offsets for each dataset, which are applied in [plot_results.R](plot_results.R). This improves the positions, but does not resolve all the location errors.

While this remains a very useful benchmark problem, these issues should be considered when interpreting the results. 

## Base results

Figure 2 shows the modelled and observed runup, plotted radially from the centre of the island. The modelled values are in orange, while the three different datasets are in green, red and black. Overall the model does a good job of predicting the observed runup heights.

Different datasets sometimes give different runup estimates (Figure 2). There are obvious location errors in some data points which plot too far offshore or inland.

The high runups near Monai (~30m) are well represented in the model (near coordinate 139.42,42.10). The test checks that runup >25 m is predicted here. With default resolution we obtain runup >26 m, which becomes closer to 30 m with further grid refinement (e.g. setting `mesh_refine=1.0_dp` in the [SWALS model](BP09.f90)). 

![Figure 2: Modelled and observed tsunami runup around Okushiri Island](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/nthmp/BP09/runup_heights_okushiri_lowresolution_omp.png)

Figure 3 shows the modelled tsunami maxima around Okushiri Island. The previous figure suggests this gives a good representation of runup onshore. There are some north-south oriented features in deeper waters associated with (spurious) jumps in the elevation data, which are not nesting related, and should be eliminated with better data.

![Figure 3: Modelled tsunami maxima near Okushiri Island](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/nthmp/BP09/max_stage_okushiri_lowresolution_omp.png)

## Consistency of the results with alternative domain partitions and local timestepping

Here we run the problem with three alternative model setups. In principle the three alternatives should lead to results that are "almost the same", if not identical. Below we check for consistency among the three approaches. They are:

1. The first run uses the domains in Figure 1 with openmp for parallelism (`openmp`).

2. The second run partitions each domain into 6 pieces as in Figure 4 (`coarray`).

    * When SWALS uses MPI (or coarray) parallelism, by default each domain is split into pieces (one per MPI rank, or coarray image). 
    * The defaults can be overridden by providing a load balance file (e.g. `md%load_balance_file = load_balance_partition.txt`). 
        * The load balance file contains one row for each of original domains. Each row includes one or more integers (one for each piece that the domain is split into). The integers give the "image index" (=`MPI rank + 1`) that runs that piece. 

3. The third run is like the second except it uses local timestepping (`omp_LocalTS`). 

    * Local timestepping allows nonlinear domains to take longer timesteps (an integer multiple of `global_dt/domain%timestepping_refinement_factor`) if that is stable according to their CFL condition.
    * No local timestepping is applied to domains using leapfrog-type solvers, as their second order accuracy relies on a fixed timestep.
    * Local timestepping is enabled by compiling SWALS with `-DLOCAL_TIMESTEP_PARTIONED_DOMAINS`, see [make_BP09_localtimestep](make_BP09_localtimestep).

![Figure 4: Domain partitioning in runs 2 and 3. Each of the original domains is split into several pieces (compare to Figure 1)](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/nthmp/BP09/elevation_okushiri_lowresolution_coarray.png)

The test code checks that runs 1, 2 and 3 give similar results. They are not bitwise identical (reasons discussed below), but the figures below show it isn't easy to notice differences in runup. For visual clarity the domain bounding boxes are not shown as partitioned in runs 2 and 3.

![Runup with Run 1](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/nthmp/BP09/runup_heights_okushiri_lowresolution_omp.png) ![Runup with Run 2](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/nthmp/BP09/runup_heights_okushiri_lowresolution_coarray.png) ![Runup with Run 3](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/nthmp/BP09/runup_heights_okushiri_lowresolution_omp_localtimestep.png)

Below we show the max-stage figures for the three runs; again we cannot notice differences. 

![Tsunami maxima with Run 1](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/nthmp/BP09/max_stage_okushiri_lowresolution_omp.png) ![Tsunami maxima with Run 2](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/nthmp/BP09/max_stage_okushiri_lowresolution_coarray.png) ![Tsunami maxima with Run 3](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/nthmp/BP09/max_stage_okushiri_lowresolution_omp_localtimestep.png)

### Time-varying model statistics in the three runs

Below we overplot time-varying model summary statistics for run 1 and run 2, and then repeat for run 1 vs run 3. The time-series are very similar in all models.

![Model statistics over time in run 1 and run 2](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/nthmp/BP09/Compare_openmp_coarray.png) ![Model statistics over time in run 1 and run 3](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/nthmp/BP09/Compare_openmp_ompLocalTS.png)

### Differences in UH in the three runs near the time of extreme runup at Monai

In the previous plots it was hard to see differences in the models, although the results were not actually identical. To better highlight differences between the models we focus on the easterly flux UH, which is more sensitive to model differences than the water surface and runup.

Figures below compare UH in runs 1 and 2 at around time=292.5s. This is near the time of Monai's extreme runup. For each domain we show both the difference plot (run 1 - run2), and the solution (run 1), using a different colour scale in each panel. While the differences are very small they are not zero. 

* More detail can be obtained by storing results in double precision rather than single precision. See the variable `output_precision` in [global_mod.f90](../../../../src/shallow_water/global_mod.f90).

![Comparison of UH in runs 1 and run 2 after 40 output time-steps (292.5s). For each domain we plot the difference, and the solution.](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/nthmp/BP09/Compare_omp_coarray_time_index_40.png)

Below we compare run 1 and run 3 in the same way. The model differences are still very small, but more noticeable than above. This is because local timestepping is more significant change to the computational method. 

![Comparison of UH in runs 1 and run 3 after 40 output time-steps (292.5s). For each domain we plot the difference, and the solution.](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/nthmp/BP09/Compare_omp_ompLocalTS_time_index_40.png)

### Differences in UH in the three runs at the end of the simulation

Over time the differences grow. To illustrate this we repeat the UH comparison (run 1 - run 2) at time=2992.5s. The differences are still small, albeit larger than at 292.5s.

![Comparison of UH in runs 1 and run 2 after 400 output time-steps (2992.5s). For each domain we plot the difference, and the solution.](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/nthmp/BP09/Compare_omp_coarray_time_index_400.png)

Below we compare run 1 and run 3 in the same way. The differences are larger, especially in nearshore domains where the model predicts eddy formation. These eddies evolve chaotically over long times, so are sensitive to details of the numerical method such as local timestepping.

![Comparison of UH in runs 1 and run 3 after 400 output time-steps (2992.5s). For each domain we plot the difference, and the solution.](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/nthmp/BP09/Compare_omp_ompLocalTS_time_index_400.png)

### Why differences are expected due to alternative domain partitions and local timestepping.

Small differences between runs 1 and 2 occur because of the domain partitioning. In finite precision floating point arithmetic, the partitioned domain coordinates are not bitwise equal to the original domain coordinates (although differences are tiny, in the least significant digits). Similarly, the flux correction algorithm leads to different summation orderings in each case. Such factors introduce tiny differences in the solutions that can grow over time. 

* It is possible to force models with different openmp/MPI configurations to give __identical__ results by prescribing the same domain partition in both cases (i.e. setting `md%load_balance_file`). This is good practice in applications. 
    * This is done in an alternative script [run_model_exact_reproduce.sh](run_model_exact_reproduce.sh). For brevity that is not run by the automated test suite, but such a test is run in [paraboloid_bowl](../../paraboloid_bowl) and in [../Hilo_Tohoku_tsunami/](../Hilo_Tohoku_tsunami).

For run 3, local timestepping changes the numerical method. Thus results differ from both run 1 and run 2. 

To further explore the sensitivity of the solution to tiny floating point differences, it is interesting to repeat run 1 while adding a tiny perturbation to the model elevations (1e-10 m). 

* This is not here, but leads to differences comparable to "run 1 - run 2". 
* The subroutine `set_initial_conditions_BP09` in [BP09.f90](BP09.f90) can be modified to do this by changing `real(dp), parameter :: random_perturbation_scale = 0.0e-10_dp` to `real(dp), parameter :: random_perturbation_scale = 1.0e-10_dp`. 

