# Analysis of model runs for scenarios ID710.5

This folder contains scripts to post-process the model runs in [../swals](../swals)

The NSW model involves multiple batches of scenarios, and herein the codes have been setup
to run scenarios ID710.5. The ID reflects the PTHA18 hazard point gauge ID. 

## Subfolders:

* [check_log_files](./check_log_files) - performs various sanity checks on model log files {mass conservation, energy decay, etc}

* [max_stage_at_a_point](./max_stage_at_a_point) - compares the Monte-Carlo solution derived from random scenarios (with a high resolution model) to the offshore PTHA solution. This is used for quality control. Broadly speaking, we expect good agreement offshore in deep water (where linear theory works well), and potentially substantial disagreement in shallow waters close to the coast.

* [probabilistic_inundation](./probabilistic_inundation) - compute inundation exceedance-rates and their epistemic uncertainties. Also compute max-stage exceedance-rates. 

* [jatwc_to_inundation](./jatwc_to_inundation) - classify scenarios according to their JATWC statistic, and thus derive inundation zones corresponding to each JATWC warning category.

* [template_data_package](./template_data_package) - documentation that is integrated into a package of output products.

## Scripts

* [make_data_package.R](make_data_package.R) is used to make a folder with output products and documentation.
