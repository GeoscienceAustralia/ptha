# Random scenario sampling of sunda2 source zone.
-------------------------------------------------

This folder contains code to sample random scenarios from the sunda2 source-zone.
Key codes are:
* [sample_random_scenarios_simple.R](sample_random_scenarios_simple.R) illustrates basic random sampling of scenarios and comparison of exceedance-rate curves with 'exact' PTHA solutions. The code is only actually used to define the non-uniform sampling effort, but it can also make nice figures.

* [select_random_scenarios.R](select_random_scenarios.R) performs stratified-importance sampling with nonuniform sampling in magnitude-bins. Both the unsegmented and segmented representations are used. The methodology follows that used for the PTHA sampling paper. 

* [generate_initial_conditions.R](generate_initial_conditions.R) makes the initial conditions for tsunami modelling.

