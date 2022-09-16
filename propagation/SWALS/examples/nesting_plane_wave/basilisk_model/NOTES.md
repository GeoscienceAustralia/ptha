The nonlinear finite-volume solvers in SWALS can be very dissipative for simple
linear plane wave propagation if we use a TVD limiter coefficient (such as
`domain%theta = 1.3`). 

I wanted to check this in another solver that uses a similar algorithm.  The
shallow water solver in Basilisk is suitable. It also uses a Kurganov type flux
function (HLLC in 2D), with stage/velocity extrapolation and limiting, and a
default limiter coefficient of 1.3.

The code [simple_wave.c](simple_wave.c) runs a plane-wave problem with Basilisk.
We find that it dissipates the plane wave and introduces assymmetries in the
solution. The behaviour is just like the SWALS FV solvers when they are
setup to use the same grid size, limiter, and a similar flux function. 

The corresponding SWALS model can be run by uncommenting code near the comment
`Alternate version` in [../nesting_reflection.f90](../nesting_reflection.f90)
(and commenting out the other definitions of those parameters).

