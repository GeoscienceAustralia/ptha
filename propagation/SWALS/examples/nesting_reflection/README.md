# Plane wave propagation with nested grids

This problem simulates 1D plane-wave propagation for a shallow water wave packet with very small amplitude (1 mm, in 100 m depth) in a periodic domain. For such a wave, linear shallow water theory should provide a very good approximation. The wave is propagated once through the domain, and the interior part of the wave packet is then compared with the initial condition. Analytically they should be the same.

The [SWALS model](nesting_reflection.f90) simulates the problem with nested grids and a range of numerical methods (staggered grid schemes `linear` and `leapfrog_nonlinear`, and finite-volume schemes `rk2` and `midpoint`). The idea is to highlight any obvious problems with particular algorithms, or the grid nesting. 

In practice we find the staggered grid algorithms provide a very accurate solution to this problem, using either the linear or the nonlinear shallow water equations. This is expected, as in general they are well suited to modelling low Froude-number flows.

![Solution with the staggered-grid `linear` flow algorithm](cycle_solution_linear.png)

![Solution with the staggered-grid `leapfrog_nonlinear` flow algorithm](cycle_solution_leapfrog_nonlinear.png)

The finite volume algorithms can also solve this problem, but are less accurate than the staggered-grid algorithms. They tend to introduce asymmetries into the numerical solution, albeit these reduce with grid refinement. 

![Solution with the finite-volume `rk2` flow algorithm](cycle_solution_rk2.png)

![Solution with the finite-volume `midpoint` flow algorithm](cycle_solution_midpoint.png)

To prevent the finite-volume slope-limiters from excessively dissipating the wave packet, this SWALS model is setup to use a non total-variation-diminishing slope-limiter (setting `md%domains(..)%theta = 4.0`; the limiters are  TVD for values less than 2). Without this, the wave packet would be noticeably attenuated. The same numerical issues arise in another model (Basilisk) that implements a very similar finite volume scheme as SWALS. 

This problem highlights that, compared to the staggered grid schemes, the finite-volume schemes are simply less well suited to long distance wave propagation at low Froude-numbers, such as global scale tsunami propagation to far-field sites. The advantages of the finite-volume schemes arise for more strongly nonlinear problems, including inundation. In practice these kinds of numerical issues can be detected using convergence tests, i.e., checking for numerical changes in the solution as the model is run with increasingly fine grid resolutions.



