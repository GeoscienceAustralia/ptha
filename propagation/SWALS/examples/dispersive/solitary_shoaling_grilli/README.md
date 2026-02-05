# Solitary wave runup with dispersion

This problem models the experiment of [Grilli et al. (1994)](http://dx.doi.org/10.1061/(ASCE)0733-950X(1994)120:6(609)) where a solitary wave runs up a linearly sloping beach.
This problem was also used in [Filippini et al. (2015)](http://dx.doi.org/10.1016/j.coastaleng.2015.02.003).

![Figure 1: Initial condition](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/dispersive/solitary_shoaling_grilli/Initial_condition.png)

The model is compared with observations at gauge 0 (which is used to temporally align the model and data) and gauges 1, 3, 5, 7, 9. 

The model reasonably reproduces the observations, with a slightly early arrival at gauges on the slope, similarly to results in Filippini et al. (2015) using Boussinesq models in flux form. 

That said, we note the initial condition used here is not a perfect solitary wave for SWALS. Furthermore, the location of gauge 9 appears to have an typo in the original paper (off by 1m), and the adjusted coordinate used here agrees with plots in Gilli et al. (1994) and Filippini et al. (2015).

## Midpoint
![Figure 2: Model vs observations, midpoint solver](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/dispersive/solitary_shoaling_grilli/Solitary_shoaling_Grilli94_midpoint.png)

## rk2
![Figure 3: Model vs observations, rk2 solver](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/dispersive/solitary_shoaling_grilli/Solitary_shoaling_Grilli94_rk2.png)

## leapfrog nonlinear
![Figure 4: Model vs observations, leapfrog nonlinear solver](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/dispersive/solitary_shoaling_grilli/Solitary_shoaling_Grilli94_leapfrog_nonlinear.png)

## Cliffs
![Figure 5: Model vs observations, cliffs solver](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/dispersive/solitary_shoaling_grilli/Solitary_shoaling_Grilli94_cliffs.png)
