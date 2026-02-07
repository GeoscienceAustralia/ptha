# NTHMP test problems 5 (dispersive version): Solitary wave on a composite beach

This problem covers benchmark 5 from the NTHMP test suite, using dispersive solvers. The test data
and a problem description is available in 
[Randy LeVeque's repository](https://github.com/rjleveque/nthmp-benchmark-problems/tree/master/BP02-DmitryN-Solitary_wave_on_composite_beach_analytic)

Three different solitary waves (cases A, B, C) are propagated over a 1D
composite beach and compared with experiments.

![Figure 1: Beach profile and gauge locations](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/dispersive/nthmp_BP05/solution_geometry_caseA_midpoint.png)

Physically, nonlinearity is significant because the wave amplitudes are not
negligible compared to the depth (especially for cases B and C). Dispersion
also matters physically because the wave length is not very long compared
to the water depth (moreso for cases B and C). Dispersion causes short waves to
form in the experiments, over regions where non-dispersive nonlinear shallow
water models predict shock formation.

For case C (and to a lesser extent case B) the dispersive terms should be
suppressed in shallower waters. The need to suppress dispersive terms in these
situations is discussed for several dispersive models in the NTHMP (2012)
workshop report. According to the latter report, the FUNWAVE model switches off
dispersive terms when the wave height exceeds 0.8 of the depth, while NEOWAVE
switches off dispersive terms when the Froude Number exceeds 0.5 (and switches
them back on when it drops below 0.15).  Herein, the model tapers the
dispersive terms to zero between threshold still-water depths, which are set
separately for each case.

The [SWALS model](BP5_testcases.f90) is setup to take the numerical method as a commandline
argument. This is used to test a few nonlinear shallow water solvers (`midpoint` and `rk2` and `leapfrog_nonlinear`).

## Nonlinear solutions - `midpoint`

Here we test the `midpoint` finite volume scheme including dispersion.

![Figure 5: Comparison of experimental and SWALS dispersive solutions for case A](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/dispersive/nthmp_BP05/solutionA_midpoint.png)

![Figure 6: Comparison of experimental and SWALS dispersive solutions for case B](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/dispersive/nthmp_BP05/solutionB_midpoint.png)

![Figure 7: Comparison of experimental and SWALS dispersive solutions for case C](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/dispersive/nthmp_BP05/solutionC_midpoint.png)

## Nonlinear solutions - `rk2`

Here we test the `rk2` finite volume scheme including dispersion. Note that when using dispersion, this may less accurate than `midpoint` because it includes the dispersive terms in a first-order accurate manner.

![Figure 5: Comparison of experimental and SWALS dispersive solutions for case A](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/dispersive/nthmp_BP05/solutionA_rk2.png)

![Figure 6: Comparison of experimental, and SWALS dispersive solutions for case B](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/dispersive/nthmp_BP05/solutionB_rk2.png)

![Figure 7: Comparison of experimental, and SWALS dispersive solutions for case C](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/dispersive/nthmp_BP05/solutionC_rk2.png)


## Nonlinear solutions - `leapfrog_nonlinear`

Here we test the `leapfrog_nonlinear` finite difference scheme including dispersion. 

![Figure 8: Comparison of experimental, and SWALS dispersive solutions for case A](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/dispersive/nthmp_BP05/solutionA_leapfrog_nonlinear.png)

![Figure 9: Comparison of experimental, and SWALS dispersive solutions for case B](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/dispersive/nthmp_BP05/solutionB_leapfrog_nonlinear.png)

![Figure 10: Comparison of experimental, and SWALS dispersive solutions for case C](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/dispersive/nthmp_BP05/solutionC_leapfrog_nonlinear.png)
