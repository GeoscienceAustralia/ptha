# Shoaling over a variable seabed

This problem was proposed by [Kennedy et al. (2000)](https://www1.udel.edu/kirby/papers/kennedy-etal-icce00.pdf) and is also
used in [Coulaud et al. (2025)](http://dx.doi.org/10.1016/j.coastaleng.2024.104645). 

A wave train is initialised with non-uniform stage, and zero velocities, over a
sloping seabed. The solution is evolved for a set time period. The initial wave periods are
only 2-3x greater than the depth, so they are very dispersive, while the wave
heights are large enough that nonlinearity also matters.

We compare the numerical water surface maxima/minima with the Whispers3D
reference solution of Coulaud et al. (2025), digitized from their paper. The
SWALS models are run at a range of resolutions to check the rate of
convergence. 

All solvers give quite similar solutions at the finest resolution. The midpoint
solver converges more quickly than the other solvers, which is not surprising
since it is the only one which retains second order accuracy for both
dispersive and nonlinear terms.

The models are also run in transposed domains (i.e. the 1D test
separately exercises the UH and VH momentum equations) to check that they 
give near identical results. The plots include both directions, but only
one will be visible since they plot almost exactly atop each other.


## Midpoint


![Figure 1: Mid point solver](FIXME)


## rk2

![Figure 2: rk2 solver](FIXME)

## Leapfrog nonlinear

![Figure 3: Leapfrog nonlinear solver](FIXME)

## Cliffs

![Figure 4: Cliffs solver](FIXME)

