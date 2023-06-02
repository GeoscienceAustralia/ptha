This folder contains code to test and post-process the random PTHA inundation scenarios. 
* [check_log_files](check_log_files) is used to sanity check model runs from their log-files.
* [max_flow_plots](max_flow_plots) was used to store some plots, and refers to related code.
* [max_stage_at_a_point](max_stage_at_a_point) is used to check consistency between the PTHA18 results and the results of our inundation simulations. The idea is that both approaches should agree fairly well in deep water sufficiently far from the coast, although this is not expected nearshore (hence why we do inundation modelling!). 
* [probabilistic_inundation](probabilistic_inundation) is used to compute exceedance-rates of inundation/max-stage/depth, including outputs that quantify epistemic uncertainties and (separately) Monte Carlo uncertainties.
* [jatwc_to_inundation](jatwc_to_inundation) is used to compute preliminary zones corresponding to JATWC warning categories. These are later manually edited (as discussed in the AJEM manuscript). 
