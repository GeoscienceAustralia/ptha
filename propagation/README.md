This folder can store codes for tsunami propagation modelling. 

Currently some Shallow-WAter-Like-SolverS (SWALS) are provided. SWALS includes a
linear-shallow-water-equations leapfrog solver which is generally suitable for
earthquake-tsunami modelling in the deep ocean, but poorly suited to modelling
nearshore flows and inundation. It also includes several
nonlinear-shallow-water-equations finite-volume solvers which are generally
suitable for modelling nearshore flows and inundation, but not well suited to
modelling global-scale propagation (i.e. inefficient compared with the
linear-leapfrog approach). SWALS supports both spherical and cartesian
coordinates; nested grids which may use different numerical schemes and
time-steps; and can be run in parallel using OpenMP and/or fortran coarrays.
Tests include a unit-test suite; a parallel unit-test suite; and a set of
tsunami model benchmark problems ( analytical / laboratory / field ).

Other well known open source shallow water solvers include GEOCLAW, JAGURS,
COMCOT, ANUGA, Basilisk, and easyWave. They vary widely in their features, and
the required computational effort.

We aim to keep the other parts of the PTHA package independent of any
particular solver. This is especially enforced in rptha, but less so for
application-specific template scripts. However most of our template scripts
could be adapted to use any solver which is able to output stage time-series at
tide gauges. In general this will require the user to change function calls
which read the solver output data, to refer to new (user provided) functions
which read from the new solver. 
