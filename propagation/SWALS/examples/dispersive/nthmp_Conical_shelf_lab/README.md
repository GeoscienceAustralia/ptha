# NTHMP currents benchmark problem 5 (dispersive version): Solitary wave over a shelf with a Conical Island

We model the propagation of a solitary wave over a shelf containing a conical island, and compare the model with measurements of the free surface and velocities at a number of gauges. Data from this problem was obtained from the [NTHMP tsunami currents benchmark suite](http://coastal.usc.edu/currents_workshop/problems.html). 

The [SWALS model](model.f90) simulates this problem on the domain shown in Figure 1. Gauges with names starting with `WG` record water-levels, while gauges `A, B` and `C` record velocities. In the experiment the wave was created using a wave-maker on the left-hand-side of the domain. Herein the [SWALS model](model.f90) simulates this by extending the domain to the left and using a solitary wave initial condition (similar to [Macais et al., 2020](https://doi.org/10.1016/j.ocemod.2020.101645)).  

![Figure 1: The problem domain and gauge locations. The wavemaker is on the left-hand-side of the domain. For this problem we emulate a wavemaker source by extending the domain to the left and using a solitary wave initial condition.](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/dispersive/nthmp_Conical_shelf_lab/domain_setup.png)

Solitary waves are defined by a balance between non-linearity and dispersion. Thus it is not surprising that the model here, including dispersion in deeper areas (tapered off between depths of 0.3 - 0.1 m), often better represents the leading waves (Figure 2) and their velocities (Figure 3) than in [the non-dispersive case](../../nthmp/Conical_shelf_lab/). This is particularly at deeper sites, where without dispersion, the leading wave tends to arrive as a shock. 

![Figure 2: Comparison of observed and modelled free-surface time-series at gauges.](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/dispersive/nthmp_Conical_shelf_lab/Stage_gauges.png)

![Figure 3: Comparison of observed and modelled velocities at gauges A, B and C. The U-velocity is in the x-direction, V-velocity in the y-direction.](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/dispersive/nthmp_Conical_shelf_lab/velocity_gauges.png)

