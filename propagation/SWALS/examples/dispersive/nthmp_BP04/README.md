# NTHMP test problem 4 (dispersive version): Solitary wave on a simple beach (experimental version)

We model the runup of a 1D solitary wave on a sloping beach and compare the model with experimental data, including multiple test cases with varying wave amplitudes and depths. The non-dispersive version of this problem is in [../../nthmp/BP04](../../nthmp/BP04). Compared to the nonlinear shallow water simulations, the dispersive model needs a slightly higher friction coefficient to accurately reproduce the runup (n = 0.008 here vs 0.005 in the other case). Dispersion is tapered to zero as the depth (below MSL) varies from 0.3 - 0.2 of the initial depth. 

This test problem is from the NTHMP benchmark suite. The test data and a problem description is available in [Randy LeVeque's repository](https://github.com/rjleveque/nthmp-benchmark-problems/tree/master/BP04-JosephZ-Single_wave_on_simple_beach). 

## The case with a small amplitude wave, $h/d = 0.0185$

The figure below compares modelled (black) and observed (red) free surface profiles for a relatively small amplitude wave. They show generally good agreement. 

![Figure 1: Comparison of modelled (midpoint) and experimental free surface at various times, low amplitude wave](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/dispersive/nthmp_BP04/Model-vs-data_0.0185.png)

## The case with a large amplitude wave, $h/d = 0.3$ 

The figure below compares modelled (black) and observed (red) free surface profiles for a relatively large amplitude wave. The dispersive model captures the shape of the wave prior to runup better than the non-dispersive model [../../nthmp/BP04](../../nthmp/BP04), because dispersion matters for this relatively large solitary wave.

![Figure 2: Comparison of modelled (midpoint) and experimental free surface at various times, high amplitude wave](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/dispersive/nthmp_BP04/Model-vs-data_0.3.png)

## Dimensionless runup vs dimensionless wave amplitude

For this case a we focus on the relationship between the runup and the initial wave amplitude, both normalised by the depth (which varies substantially among the experiments). With our chosen Manning coefficient of 0.008, the `midpoint` solver shows good agreement with the general trend of the experimental data. The Manning coefficient was tuned to give reasonably good results, and is different to the value (0.005) that gave good results for the non-dispersive model (albeit both are characteristic of low-friction surfaces).

![Figure 3: Dimensionless runup as a function of dimensionless wave amplitude in experiments and model (midpoint).](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/dispersive/nthmp_BP04/Runup_scaling_plot.png)
