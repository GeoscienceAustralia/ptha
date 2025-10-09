# NTHMP test problem 6 (dispersive version): Solitary wave runup on a conical island.

We model the runup of three different solitary waves around a conical island, including dispersion in the solver (non-dispersive version [here](../../nthmp/BP06)). The three cases (A, B, and C) are forced with increasingly large solitary wave initial conditions, slightly less than 5%, 10% and 20% of the offshore depth. They were studied experimentally by [Briggs et al., 1995](https://doi.org/10.1007/bf00874384) and have been modelled in many publications.

The test problem is from the NTHMP benchmark suite. The test data and a problem description is available in [Randy LeVeque's repository](https://github.com/rjleveque/nthmp-benchmark-problems/tree/master/BP06-FrankG-Solitary_wave_on_a_conical_island). This includes gauge time-series between the wavemaker and the island (gauges 1,2,3,4), gauge time-series around the island (gauges 6, 9, 16, 22), and records of the runup maxima around the island. Note the locations of gauges 1-4 differ for cases A, B, and C.

![Conical Island flume geometry and gauge locations. The wavemaker is near the left edge of the domain. Gauge locations 1-4 vary slightly between cases A, B and C](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/dispersive/nthmp_BP06/Flume_plot_A_default.png)

The [SWALS model](BP06.f90) uses a nested grid around the island. The wave forcing was derived from an analytical initial condition for a solitary wave. Dispersion is used in deeper waters, and linearly tapered to zero at depths (below MSL) between 0.1 - 0.05 m.

## Differing runup observations are reported for this problem

For this problem the runup data provided for the NTHMP benchmark, and used in many publications (e.g.  [Tolkova, 2014](https://doi.org/10.1007/s00024-014-0825-8) and [Horillo et al., 2014](10.1007/s00024-014-0891-y)), differs noticeably to data in other publications ([Liu et al., 1995](https://doi.org/10.1017/S0022112095004095), [Choi et al., 2007](https://doi.org/10.1016/j.coastaleng.2007.02.001), [Ma et al., 2019](https://doi.org/10.1080/19942060.2019.1642960)). The latter publications also attribute their data to [Briggs et al., 1995](https://doi.org/10.1007/bf00874384), and there do not appear to be any other differences in the problem description. Runup maxima are especially different for Case A, where runups in the latter publications are consistently 30-50% larger than for the NTHMP data.

Comparisons with models and data in the literature do not consistently show better agreement with one or other dataset, but there is considerable variation between models. The reasons are not clear, but several factors likely contribute.

* The wave forcing in this problem can be defined several ways. Some studies use empirical wavemaker data. Others specify an initial solitary wave (as for FUNWAVE herein). Others adapt the measured wave time-series to specify a boundary condition (as for SWALS herein).

* The experiment involves solitary waves, which feature both dispersion and nonlinearity. They are not well represented with non-dispersive shallow water models (used in many studies). However some shallow water models can emulate physical dispersion by deliberately using a coarse grid (e.g. [Tolkova, 2014](https://doi.org/10.1007/s00024-014-0825-8)). 

But even when this problem is approached with more complicated physics (i.e. dispersive or 3D models), we see reports of good agreement for some models using the NTHMP dataset (e.g. SELFE in the [NTHMP 2011 workshop report](https://nws.weather.gov/nthmp/documents/nthmpWorkshopProcMerged.pdf)), and for other models using the alternative runup dataset (e.g. [Choi et al., 2007](https://doi.org/10.1016/j.coastaleng.2007.02.001), [Ma et al., 2019](https://doi.org/10.1080/19942060.2019.1642960)). Because it is unclear which dataset should be used, both are included below. 

## Reference FUNWAVE model, without dispersion, using analytical forcing

We also compare the SWALS runup results with a high-resolution FUNWAVE simulation. FUNWAVE was run __without dispersion__ on a 2.5cm uniform grid, and initialised with an analytical solitary wave. [See here for FUNWAVE model setup files](funwave_comparison). 

This serves as a high-order reference solution to the nondispersive shallow water equations. It is not intended to test the FUNWAVE model itself, for which better results might be obtained using dispersion.

# Results for Case A

Figure 1 compares the SWALS model with the gauge time-series. For early times the model agrees well with the offshore gauges (1-4) and the gauges around the island (6,9,16,22). At later times the waves are affected by interaction with the flume boundary, which is designed to reduce reflections, and not perfectly represented in the model (which here uses reflective north/south boundaries and a Flather east/west boundary). 

![Figure 1: Modelled and observed time-series at offshore gauges, Case A](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/dispersive/nthmp_BP06/Gauges_plot_A_default.png)

Figure 2 compares the modelled and observed runups around the island. Note the large difference between the two runup datasets. With dispersion, the SWALS runups are slightly smaller than the reference FUNWAVE shallow water simulation and the "alternative" runup data. Both modelled runups are significantly larger than the NTHMP runup data.

* If dispersion is included in the FUNWAVE simulation then the modelled runup reduces also slightly (in between the two datasets).

![Figure 2: Modelled and observed runup maxima around the island, Case A](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/dispersive/nthmp_BP06/Runup_plot_A_default.png)

# Results for Case B

Figure 3 compares the SWALS model with the gauge time-series. 

![Figure 3: Modelled and observed time-series at offshore gauges, Case B](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/dispersive/nthmp_BP06/Gauges_plot_B_default.png)

Figure 4 compares the modelled and observed runups around the island. 

![Figure 4: Modelled and observed runup maxima around the island, Case B](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/dispersive/nthmp_BP06/Runup_plot_B_default.png)

# Results for Case C

Figure 5 compares the SWALS model with the gauge time-series.  

![Figure 5: Modelled and observed time-series at offshore gauges, Case C](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/dispersive/nthmp_BP06/Gauges_plot_C_default.png)

Figure 6 compares the modelled and observed runups around the island.

![Figure 6: Modelled and observed runup maxima around the island, Case C](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/dispersive/nthmp_BP06/Runup_plot_C_default.png)
