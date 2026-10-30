# Post processing of sampled scenarios.

Key directories are:
* [gauge_plots](gauge_plots) validates the model with a source inversion against tide gauge records of historic earthquake tsunamis.
* [max_stage_at_a_point](max_stage_at_a_point) compares the Monte Carlo results to PTHA18 at offshore sites, to check it all worked.
* [probabilistic_inundation](probabilistic_inundation) computes inundation rates, and thresholds of flow variables with given exceedance rates. Ensure [max_stage_at_a_point](max_stage_at_a_point) is done first.
* [jatwc_to_inundation](jatwc_to_inundation) defines model-based zones associated with each of the JATWC warning zones.
* [convert_max_stage_outputs_to_AHD](convert_max_stage_outputs_to_AHD) translates a few max-stage flow variables into an AHD vertical datum (versus the original datum of "m above the local high tide from TPXO9v5a").
