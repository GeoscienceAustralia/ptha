Compare the PTHA18 maximum-stage at a given offshore site with the result of the nonlinear SWALS models. This is a useful way to check the high-resolution Monte Carlo results. Basically the latter should agree well with the offshore PTHA in deep water far from the coast (i.e. at sites where limitations in the offshore PTHA's hydrodynamic model are not important). But at sites closer to the coast or in shallower water, we do not expect the offshore PTHA's hydrodynamic model to work so well -- and we expect differences.

The key code is [extract_max_stage_at_a_point_extended.R](extract_max_stage_at_a_point_extended.R)
* The script [run_a_few.sh](run_a_few.sh) applies it to a number of cases. 

