# Modelled inundation products corresponding to JATWC coastal zones and threat levels.

The Joint Australian Tsunami Warning Centre (JATWC) provides tsunami warnings,
with separate warnings provided in each coastal zone. As of November 2024 there
are three possible warning categories (threat levels):
* No threat
* Marine warning
* Land warning


## Contents

* `ATWS_COASTAL_ZONE_POLYGONS/` contains a shapefile with the ATWS coastal zones, each of which will receive its own warning category.
  * This shapefile was kindly provided by Robert Greenwood (Bureau expert on JATWC tsunami warning system).
  * Coastal zone names in the `ATWS_ZONES` attribute correspond to zone names used for the output products.

* `Inundation_zones` contains one sub-folder for each ATWS coastal zone in NSW including Lord Howe Island (and Norfolk Island). For each coastal zone, the products are:
  * Polygons depicting the model-based inundation zones for `no threat`, `marine warning` and `land warning` categories. 
    * These are in folders with names like `Sydney-Coast_land_warning_with-PTHA-exrate-limit_84pc_4e-04`
      * For this example the coastal zone is `Sydney Coast` and the warning category is `land warning`.
    * All zones are limited to areas with non-zero inundation at an exceedance-rate of 0.0004 (=1/2500) with 84th percentile epistemic uncertainty.
      * This is a rare event, even assuming a conservative earthquake frequency model.
      * It serves as a probabilistic alternative to a "credible worst case scenario" (the latter concept is otherwise difficult to define for tsunamis).
    * The output polygons extend slightly beyond their coastal zone. 
      * This is useful to provide insights on possible effects in nearby zones (e.g. when considering the boundaries between zones). 
      * But it also means that, if neighbouring coastal zone polygons are displayed together, their polygons may have some overlap. 
      * If this becomes a problem then we can do further processing to make the boundaries sharp.
  * Rasters depicting the `maximum stage` (waterlevel in m AHD) over the combination of `no_threat` scenarios,  and (separately) `marine_warning` scenarios 
    * For instance the folder `Sydney-Coast_marine_warning_max_stage_tifs` contains the maximum stage (waterlevel) reached by ANY of our marine warning scenarios in the `Sydney-Coast` zone.
    * These results are only provided for `no_threat` and `marine_warning`. 
      * There is no similar result for `land_warning` scenarios, because they should be limited to a the maximum stage corresponding to a rate of 1/2500 with 84th percentile epistemic uncertainty. 
      * Instead we suggest you look for the maximum stage at 1/2500 84th percentile within in the `../Probabilistic_inundation_hazard/` folder.

