# 1D dam-break problem

This problem simulates a 1D dam-break over uniform topography using the frictionless nonlinear shallow water equations. The analytical solution to this problem is well known, for example [Wu et al., 1999](https://ascelibrary.org/doi/10.1061/%28ASCE%290733-9429%281999%29125%3A11%281210%29).

The [SWALS model](dam_break.f90) is setup to take the upstream and downstream water depths as command-line arguments. In this test the upstream depth is always 1m, while the downstream depth takes four different values (0.5m, 0.1m, 0.01m, 0.0001m). The modelled stage and velocity at a time of 15s after the dam-break are compared with the analytical solution, and show good agreement in all cases. 

![Initial downstream depth of 0.5m](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/dambreak/dam_break_numerical_vs_analytical_H0_0.5_H1_1.png)

![Initial downstream depth of 0.1m](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/dambreak/dam_break_numerical_vs_analytical_H0_0.1_H1_1.png)

![Initial downstream depth of 0.01m](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/dambreak/dam_break_numerical_vs_analytical_H0_0.01_H1_1.png)

![Initial downstream depth of 0.0001m](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/dambreak/dam_break_numerical_vs_analytical_H0_1e-04_H1_1.png)
