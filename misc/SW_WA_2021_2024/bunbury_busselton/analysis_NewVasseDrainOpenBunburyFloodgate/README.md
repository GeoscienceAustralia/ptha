# Analyses of the random PTHA scenarios

This folder is similar to contents in ../analysis, but the codes are more updated, and they were used to process a model variant with an open floodgate at Bunbury and improved elevation data around the Vasse Diversion drain

This folder contains:

* `check_log_files` - performs various sanity checks on model log files {mass conservation, energy decay, etc}

* `max_stage_at_a_point` - compares the Monte-Carlo solution derived from random scenarios (with a high resolution model) to the offshore PTHA solution. Broadly speaking, we expect good agreement offshore in deep water (where linear theory works well), and potentially substantial disagreement in shallow waters close to the coast.

* `probabilistic_inundation` - compute inundation exceedance-rates and their epistemic uncertainties. Also compute max-stage exceedance-rates. 

* `jatwc_to_inundation` - classify scenarios according to their JATWC statistic, and thus derive inundation zones corresponding to each JATWC warning category.