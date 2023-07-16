Compare the PTHA18 maximum-stage at a given offshore site with the result of the nonlinear SWALS models. This is a useful way to check the high-resolution Monte Carlo results. Basically the latter should agree well with the offshore PTHA in deep water far from the coast (i.e. at sites where limitations in the offshore PTHA's hydrodynamic model are not important). But at sites closer to the coast or in shallower water, we do not expect the offshore PTHA's hydrodynamic model to work so well -- and we expect differences.

The key code is [extract_max_stage_at_a_point.R](extract_max_stage_at_a_point.R). The script [run_a_few.sh](run_a_few.sh) applies it to a number of cases. 

Separately there is a code to check a particular exceedance-rate at Hillarys ([this one](check_hillarys.R)), and some other old plotting code to [compare offshore and onshore stage values](plot_offshore_vs_onshore_stage.R).
