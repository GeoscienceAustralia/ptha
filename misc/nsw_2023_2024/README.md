# Files for NSW broadscale model

These files were used for some NSW wide modelling in 2023-2024.

Many were derived by modifying files in [../SW_WA_2021_2024/greater_perth_revised2023/](../SW_WA_2021_2024/greater_perth_revised2023/). For this project we use a different scenario sampling method (importance sampling, rather than stratified importance sampling) and so many of the scenario creation and probabilistic calculation scripts have changed.  

The folders are:
* `analysis_multiple_importance_sampling` - Final calculations that combine results from each batch of scenarios 
* `analysis_scenarios_ID4186.3` - Probabilistic calculations for the third batch of scenarios, which can be run after the models in `swals` have been run.
* `analysis_scenarios_ID1315.5` - As above for the second batch of scenarios
* `analysis_scenarios_ID710.5` - As aboev for the first batch of scenarios
* `breakwalls` - Input data on breakwalls that are burned into the swals model elevation
* `elevation` - Information on elevation data used for the swals model
* `inverts` - Input data on inverts that are burned into the swals model elevation.
* `multidomain_design` - Information on the multidomain structure of the swals model
* `point_gauges` - locations where the model stores gauge outputs
* `sources` - Tsunami initial conditions for historic and hazard scenarios
* `swals` - Main hydrodynamic calculations.
