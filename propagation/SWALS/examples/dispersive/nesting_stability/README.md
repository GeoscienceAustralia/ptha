# Nesting (in)stability when using dispersion

Stationary flow in deep water with irregular topography and a nested
grid. The grid size is significantly finer than the typical depth. 

Analytically the solution should be stable for all time, however, round-off
level flow will occur due to the nested finite-volume domain. Nesting
algorithms can amplify this over time, which can be particularly severe when
using dispersion. With enough time, problematic cases develop significant
speeds.

Thus, this model is useful for detecting unstable nesting algorithms. 

We test cases where the outer grid uses a staggered solver, and a cell centred
solver, with and without MPI (which causes additional partitioning of the
domains). The inner grid always uses the cell centred midpoint solver.

By default the test is not run for long enough to strongly demonstrate the
instabilities (to keep the run time short enough). Instead we use a PASS/FAIL
criterion based on the maximum speed, which can distinguish some problematic
models that we tested. 

If trying out other algorithms I'd suggest running the models for longer (~10x)
and interactively checking the results. 
