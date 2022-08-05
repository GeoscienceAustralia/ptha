# NTHMP test problem 6: Solitary wave runup on a conical island.

We model the runup of three different solitary waves around a conical island. The three cases (A, B, and C) are forced with increasingly large solitary wave amplitudes, slightly less than 5%, 10% and 20% of the offshore depth. They were studied experimentally by [Briggs et al., 1995](https://doi.org/10.1007/bf00874384) and have been modelled in many publications.

The test problem is from the NTHMP benchmark suite. The test data and a problem description is available in [Randy LeVeque's repository](https://github.com/rjleveque/nthmp-benchmark-problems/tree/master/BP06-FrankG-Solitary_wave_on_a_conical_island). This includes gauge time-series near the island, records of the runup maxima around the island, and details of the gauge locations (which vary between cases A, B, and C).

The [SWALS model](BP06.f90) simulates this problem using a nested grid. The wave forcing is created by modifying the observed wave time-series (which are interior to the model domain) to estimate how the boundary should be forced; the approach is heuristic but shown in [convert_obs_to_wavemaker.R](convert_obs_to_wavemaker.R). The default nonlinear solver is used (`rk2`) by passing `default` as the second command line argument. The first command line argument takes the value 1, 2 or 3, corresponding to cases A, B and C.

The [SWALS model code](BP06.f90) has a variety of options that allow use of other solvers, grid resolutions and forcings, reflecting experiments with this problem. For example the solitary wave forcing can be created with an analytical initial condition; an approximate wavemaker forcing; or (as used herein) a forcing based on the observed wave time-series. A model setup using the `cliffs` solver was implemented to mimic the results of [Tolkova, 2014](https://doi.org/10.1007/s00024-014-0825-8). While these variations are not discussed in detail here, they are useful for exploring the sensitivity of results to details of the model setup. 

# Note on the runup data
For this problem the runup data provided for the NTHMP benchmark, and used in many publications (e.g.  [Tolkova, 2014](https://doi.org/10.1007/s00024-014-0825-8) and [Horillo et al., 2014](10.1007/s00024-014-0891-y)), differs noticeably to data in other publications: ([Liu et al., 1995](https://doi.org/10.1017/S0022112095004095), [Choi et al., 2007](https://doi.org/10.1016/j.coastaleng.2007.02.001), [Ma et al., 2019](https://doi.org/10.1080/19942060.2019.1642960)). The latter publications also attribute their data to [Briggs et al., 1995](https://doi.org/10.1007/bf00874384), and there do not appear to be any other differences in the problem description. Runup maxima are especially different for Case A, where runups in the latter publications are consistently 30-50% larger than for the NTHMP data.

Comparisons with models and data in the literature do not consistently show better agreement with one or other dataset, but there is considerable variation between models. The reasons are not clear, but several factors likely contribute.

One factor is that alternative approaches are used to define the wave forcing in this problem. Some studies use empirical wavemaker data. Others specify an initial solitary wave (as herein). Others adapt the measured wave time-series to specify a boundary condition.

Another factor is that solitary waves feature both dispersion and nonlinearity. Thus they are not well represented with shallow water models (used in many studies). Some shallow water models can emulate the physical dispersion with numerical dispersion by deliberately using a coarse grid (e.g. [Tolkova, 2014](https://doi.org/10.1007/s00024-014-0825-8)). Yet even when this problem is approached with more complicated physics (i.e. dispersive or 3D models), we see reports of good agreement in case A for some models using the NTHMP dataset (e.g. SELFE in the [NTHMP 2011 workshop report](https://nws.weather.gov/nthmp/documents/nthmpWorkshopProcMerged.pdf)), and for other models using the alternative runup dataset (e.g. [Choi et al., 2007](https://doi.org/10.1016/j.coastaleng.2007.02.001), [Ma et al., 2019](https://doi.org/10.1080/19942060.2019.1642960)).

Because it isn't clear to the author which dataset should be used, herein we compare the model with both sets of data. We also compare to a high-resolution simulation using a reference FUNWAVE model, which was run without dispersion on a 2.5cm uniform grid, and initialised with an analytical solitary wave. ( [See here for FUNWAVE model setup files](funwave_comparison) ). This serves as a high-order reference solution to the shallow water equations, while using a different style of wave forcing than the SWALS model. 

# Results for Case A

Figures below compare the SWALS model with the gauge time-series, and the runup around the island. The model agrees well with the offshore gauges at early times. At later times the waves are affected by interaction with the flume boundary, which is designed to reduce reflections, and not perfectly represented in the model (which uses a Flather boundary so waves can radiate from the domain). The runups agree well with the reference FUNWAVE shallow water simulation, and with the "alternative" runup data, but are significantly larger than the NTHMP runup data.

![Figure 1: Modelled and observed time-series at offshore gauges, Case A](Gauges_plot_A_default.png)

![Figure 2: Modelled and observed runup maxima around the island, Case A](Runup_plot_A_default.png)

# Results for Case B

Figures below compare the SWALS model with the gauge time-series, and the runup around the island. The model shows reasonable agreement with the offshore gauges. The runups are similar to the reference FUNWAVE shallow water simulation, but less so than in case A because of the different model forcing. If SWALS uses the analytical forcing then it better agrees with the FUNWAVE result. The SWALS model also agrees reasonably well with the data. 

![Figure 3: Modelled and observed time-series at offshore gauges, Case B](Gauges_plot_B_default.png)

![Figure 4: Modelled and observed runup maxima around the island, Case B](Runup_plot_B_default.png)

# Results for Case C

Figures below compare the SWALS model with the gauge time-series, and the runup around the island. The runups agree quite well with the reference FUNWAVE shallow water simulation and both datasets. 

![Figure 5: Modelled and observed time-series at offshore gauges, Case C](Gauges_plot_C_default.png)

![Figure 6: Modelled and observed runup maxima around the island, Case C](Runup_plot_C_default.png)
