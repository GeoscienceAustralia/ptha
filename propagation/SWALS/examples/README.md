Test problems and example model setups
--------------------------------------

You can run all the test problems in one hit using the script in [../tests/validation_tests/][../tests/validation_tests/], which takes about 35 minutes on the authors home desktop. Before trying this, make sure you can get the unit-tests and parallel-tests to PASS (as discussed [here](../README.md)). In each subdirectory this will compile and run the model, make a series of the plots, and report multiple PASS/FAIL criteria for the problem. 

All tests should PASS for the default setup. However if you are trying a non-default flow algorithm it is possible for some of the tests to FAIL even if the results are ok (because some of the tests require ad-hoc threshold criteria, which can arguably be relaxed). In that case examination of the plots will give insight into the model performance.

The problems are:

* circular_island -- This compares the model with a 2D analytical solution for the linear shallow water equations, for periodic waves around a conical island.

* dambreak -- This compares the model with analytical solutions to the dam-break problem with various initial depths.

* generic_example -- This is useful for modelling single-grid global-scale tsunami propagation. It is also set up as a regression test; we simulate a PTHA scenario that is loosely similar to the Tohoku tsunami, compare a range of flow algorithms with observations, and check that the results are the same as when we first set up the problem (to within some tolerance that will likely catch changes in the model behaviour, but should allow for differences between compilers). We also check that adding a uniform rise-time to the earthquake source is equivalent to smoothing time (as is mathematically true for linear equations).

* isolated_building -- This compares the model with a well-know dam-break type experiment.

* landslide_tsunami -- This compares the model with a 1D analytical solution to runup over a plane beach, with an initial condition that is heuristically like a landslide.

* merewether -- This compares the model with a flooding scenario for merewether, similar to [the ANUGA test](https://github.com/GeoscienceAustralia/anuga_core/tree/master/validation_tests/case_studies/merewether), which is based on field observations and laboratory experiments.

* nesting_reflection -- This is a simple 1D plane wave propagation test problem which includes a nested grid. 

* parabolic_canal -- This compares the model with a 2D analytical solution for flow in a parabolic canal, which also includes a coriolis term and periodic wetting and drying.

* paraboloid_bowl -- This compares the model with a 2D analytical solution for flow in a paraboloid bowl, which is periodic with strong wetting and drying.

* periodic_convergence -- This checks the order of accuracy of the model via convergence tests, using test problem from the literature.

* periodic_multidomain -- This is not a test, but illustrates a multidomain model with periodic east-west boundary conditions, and allows a rise-time in the earthquake source.

* radial_dam_break -- This runs a radial dam break problem.

* uniform_channel -- This tests the model with steady-uniform-flow in a uniform channel

* uniform_slope -- This tests the model with steady-uniform-flow in a uniform slope

* uniform_slope_shallow -- This tests the model with steady-uniform-flow in a uniform slope, when the flow is very shallow.

See [here](./nthmp) for additional tests using the NTHMP benchmark problems, which are also run by the aforementioned validation test script.
