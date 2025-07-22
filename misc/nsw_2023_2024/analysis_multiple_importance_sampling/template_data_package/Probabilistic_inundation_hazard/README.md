# Modelled tsunami inundation hazard products

This folder contains model outputs related to:
* The long-term average frequency of tsunami inundation (events/year)
* Tsunami intensity metrics, such as the maximum depth, or the maximum stage (water-level) that are exceeded with a given long-term average frequency.
* Flood hazard categories from the [NSW flood risk management guideline FB03](https://www.environment.nsw.gov.au/research-and-publications/publications-search/flood-hazard)

The long-term average frequencies of large earthquakes are uncertain due to
limitations of current day scientific knowledge (termed "epistemic
uncertainty"). As a result the tsunami frequencies are also subject to
epistemic uncertainty. Herein the epistemic uncertainty can be described in
various ways, including using percentiles (e.g. the 84th percentile of the
epistemic uncertainty distribution) or an average.

All results have been derived by combining the 2018 Australian Probabilistic
Tsunami Hazard Assessment with inundation models that employ some conservative assumptions:
* The tsunami inundation is modelled without tides, assuming a constant background sea-level of 1.1m AHD (which is similar to highest astronomical tide for much of the NSW open coast).
* Constant friction (Manning's n of 0.03) on a bare-earth DEM. This is reasonable for relatively open ground, but does not account for the additional drag expected from buildings or dense vegetation.

This conservatism suggests our modelled tsunami inundation frequencies and
depths will be "on the high side" and should be kept in mind when interpreting
the results. For emergency management applications this is typically a good
thing.


## Contents

* `exceedance_rate_of_inundation_1mm_depth` contains rasters depicting the modelled inundation frequency (events/year).
  * These products illustrate how often inundation might occur.
  * If the modelled depth exceeds 1mm, then we say inundation has occurred. 
  * Multiple products are provided to represent epistemic uncertainties in the inundation frequency
    * `logic_tree_mean` contains results for the average earthquake frequency model
    * `84th_percentile` contains results for the 84th percentile earthquake frequency model (more conservative, higher hazard)
    * `16th_percentile` contains results for the 84th percentile earthquake frequency model (less conservative, lower hazard)
  * For each case we have a set of rasters (tifs)
    * There is 1 raster per domain in `../domains_shapefile` 
      * EXCEPTION: We exclude coarse-resolution domains in some cases because they are not needed for applications and we can avoid some heavy calculations.
      * Inundation zones are always provided
    * The rasters may overlap, but at any point only 1 raster will contain non-missing data.
    * It is convenient to treat these files as 1 using a VRT file.
      * This can be created, e.g. with `gdalbuildvrt -resolution highest name_of_vrt_file_to_create.vrt string_matching_tif_*_files.tif`
      * I'm not sure whether ESRI products have good support for VRT. If not, it may be better to build a "Raster Mosaic" (?).

* `flow_variables_1in2500_at_84th_percentile` contains flow variables with an exceedance-rate of 0.0004 (=1/2500) at the 84th percentile epistemic uncertainty. These illustrate the tsunami intensity for "a rare event according to a conservative earthquake frequency model".
  * `depth_above_initial_condition_with_exceedance_rate_1in2500_at_84th_percentile_masked_below_MSL` shows the modelled depth above the initial condition. This is the same as the depth for sites with elevation above 1.1m, while for sites with lower elevation, it is the maximum water elevation minus 1.1m. Sites with elevation below MSL have been removed.
    * We report depth-above-initial-condition to avoid exaggeration of tsunami impacts at low-lying sites (that were flooded by the initial condition).
    * Depth values are calculated to within a few mm tolerance. 
    * Sites with `elevation < zero` or `depth < 1.5mm above initial condition` have been masked for ease of interpretation.
    * We only provide results for high-resolution domains to reduce the computational load. 
      * This includes all NSW inundation sites.
  * `max_stage_with_exceedance_rate_1in2500_at_84th_percentile` 
    * This shows the tsunami maximum-stage (i.e. maximum water-level in m AHD)
    * It is clipped below 1.1015m (very slightly above 1.1m to prevent ambiguous interpretation of dry areas) 
    * To ease interpretation we mask sites that are dry, or where the tsunami is extremely small (much less than 1 mm).
    * We only provide results for high-resolution domains to reduce the computational load. 
      * This includes all NSW inundation sites.
  * `max_flux_with_exceedance_rate_1in2500_at_84th_percentile` 
    * This shows the tsunami maximum-flux (i.e. speed x depth, with units of in m^2/s)
    * It is clipped below 0.015m2/s to prevent ambiguous interpretation of dry areas
    * To ease interpretation we mask sites that are dry, or where the tsunami is extremely small (much less than 1 mm).
    * We only provide results for high-resolution domains to reduce the computational load. 
      * This includes all NSW inundation sites.
  * `max_speed_with_exceedance_rate_1in2500_at_84th_percentile` 
    * This shows the tsunami maximum-speed (i.e. with units of in m/s)
    * It is clipped below 0.0015 m/s to prevent ambiguous interpretation of dry areas
    * To ease interpretation we mask sites that are dry, or where the tsunami is extremely small (much less than 1 mm).
    * We only provide results for high-resolution domains to reduce the computational load. 
      * This includes all NSW inundation sites.
  * `flood_hazard_categories_1in2500_at_84th_percentile`
    * This shows the flood hazard categories from the NSW flood risk management guideline FB03
    * It is was computed using the max-flux variable, with limits imposed by max-speed and depth-above-initial-condition. The latter was used to downweight tsunami impacts at sites that are spuriously flooded by the initial sea level. 
    * We only provide results for high-resolution domains to reduce the computational load. 
      * This includes all NSW inundation sites.

* `flow_variables_1in250_at_50th_percentile` contains flow variables with an exceedance-rate of 0.004 (=1/250) at the 50th percentile epistemic uncertainty. These illustrate the tsunami intensity for "a non-catastrophic event, generally larger than the Chile 1960 tsunami".
  * The outputs are completely analogous to those described above (for 1/2500 at 84th percentile), but for 1/250 at 50th percentile.


