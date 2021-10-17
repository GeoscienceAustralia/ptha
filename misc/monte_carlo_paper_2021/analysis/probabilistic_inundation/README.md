# Analysis of Monte-Carlo inundation model results
------------------------------------------------

The codes here can be run after all simulations in the [../../swals](../../swals) folder have been run. They do various probabilistic calculations with the high-resolution model results, and make Figures 10, 11, 12 of the paper. 

## Key codes to run the main computations and make Figures 10 and 11:

* [run_make_probabilistic_inundation.sh](run_make_probabilistic_inundation.sh) submits the script [make_probabilistic_inundation.sh](make_probabilistic_inundation.sh). This calls a number of the R scripts that compute most of the onshore hazard results (everything except rasters depicting the 16/84th percentile epistemic uncertainties).

* [run_make_depth_epistemic_uncertainty_rasters.sh](run_make_depth_epistemic_uncertainty_rasters.sh) computes rasters depicting the 2%-in-50-years 84th percentile depths.

* [run_make_depth_epistemic_uncertainty_rasters_msl0.8.sh](run_make_depth_epistemic_uncertainty_rasters_msl0.8.sh) computes rasters depicting the 2%-in-50-years 84th percentile depths assuming the tsunami maxima coincides with a sea-level of 80cm (at the time of writing this is a typical monthly tidal maxima at Nuku'alofa).

* [run_make_depth_epistemic_uncertainty_rasters_16pc.sh](run_make_depth_epistemic_uncertainty_rasters_16pc.sh) computes rasters depicting the 2%-in-50-years 16th percentile depths.


## Code used to make Figure 12 

* [create_base_elevation_rasters.R](create_base_elevation_rasters.R) makes raster files from the model outputs, that are useful for plotting. 
* [create_openstreetmap_land_water_rasters.R](create_openstreetmap_land_water_rasters.R) makes rasters that are used to clip the inundation hazard results to inside the openstreetmap coastline.
* [plot_depth_1_in_2475.R](plot_depth_1_in_2475.R) makes Figure 12 in the paper, assuming the previous codes have been run.

## Other 

* [dominance_of_kermadectonga2_source.R](dominance_of_kermadectonga2_source.R) computes the fraction of exceedance-rates at Site P that come from the kermadectonga2 source-zone.
