# Multiple importance sampling calculations

This folder combines the calculations for three separate batches of scenarios
using multiple importance sampling (MIS). It includes the following folders which
mirror those in `../analysis_scenarios_IDxxxx.x`, but with different details to
enable MIS.
* `max_stage_at_a_point` compares the MIS calculations with the offshore PTHA
* `probabilistic_inundation` computes the MIS logic-tree-mean inundation rate, 84th/16th percentile inundation rates, various flow statistics at thresholds, and a flood hazard categorisation.
* `jatwc_to_inundation` computes zones corresponding to JATWC tsunami warning categories  

In addition there is a script `make_data_package.R` which combines specific outputs with documentation in `template_data_package` to create a folder with the final outputs.

There is also a folder `extra_plots` with code used to make some figures for the report, which might be useful to adapt in other projects. 
