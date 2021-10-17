Analysis of Monte-Carlo inundation model results
------------------------------------------------

The codes here can be run after all simulations in the [../../swals](../../swals) folder have been run.

Key codes are:

* [run_make_probabilistic_inundation.sh](run_make_probabilistic_inundation.sh) submits the script [make_probabilistic_inundation.sh](make_probabilistic_inundation.sh). This calls a number of the R scripts that compute onshore hazard results.

* [make_depth_epistemic_uncertainty_rasters.sh](make_depth_epistemic_uncertainty_rasters.sh) computes rasters depicting the 2%-in-50-years 84th percentile depths.

* [make_depth_epistemic_uncertainty_rasters_msl0.8.sh](make_depth_epistemic_uncertainty_rasters_msl0.8.sh) computes rasters depicting the 2%-in-50-years 84th percentile depths assuming the tsunami maxima coincides with a sea-level of 80cm (at the time of writing this is a typical monthly tidal maxima at Nuku'alofa).

* [make_depth_epistemic_uncertainty_rasters_16pc.sh](make_depth_epistemic_uncertainty_rasters_16pc.sh) computes rasters depicting the 2%-in-50-years 16th percentile depths.
