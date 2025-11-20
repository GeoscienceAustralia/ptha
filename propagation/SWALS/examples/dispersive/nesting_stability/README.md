# Nesting (in)stability when using dispersion

A model of stationary flow in deep water with irregular topography and a nested
grid, where the grid size is significantly finer than the typical depth. The nested
domain uses the midpoint solver, while the outer domain solver type can be set on 
the command line.

Analytically the solution should be stable for all time, however, round-off
level flow will occur due to the nested finite-volume domain. Nesting
algorithms can amplify this over time, which can be particularly severe when
using dispersion. Thus, the model is useful for detecting unstable nesting
algorithms.

The test is not run for long enough to strongly demonstrate the instabilities.
Instead we use a PASS/FAIL criterion based on the maximum speed, which can
distinguish some problematic models that we tested. 

If trying out other algorithms I'd suggest running the models for longer (~10x)
and interactively checking the results. With enough time, problematic cases
will develop significant speeds.
