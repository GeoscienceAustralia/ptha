# Compute inundation zones corresponding to JATWC warning categories, with upper limits based on PTHA.

After the scripts have been run, the key outputs are polygon shapefiles in folders named like `Inundation_zones/{zone_name}/{zone_name}_{warning_type}_with_PTHA-exrate-limit_{...}/`

# Details of folder contents (once the scripts have been run)

This folder contains the following data, along with some related scripts:

* ./ATWS_ZONES 
    * Shapefile with the ATWS warning zones, provided by Robert Greenwood from BOM (who is a technical expert in the JATWC). 
    * During a tsunami warning, each zone is assigned a separate threat level by JATWC. 

* ./Inundation_zones 
    * Onshore inundation zones for different threat levels in different ATWS warning zones, computed from our inundation scenarios, limited by a very rare PTHA exceedance-rate.
    * The key outputs are polygon shapefiles, in folders named like 
        * `Inundation_zones/{zone_name}/{zone_name}_{warning_type}_with_PTHA-exrate-limit_{...}/`
    * For each modelled scenario, in each warning zone, the threat level is defined with a statistic H (representing the 95th percentile of the tsunami maximum stage in the zone, where depths > 20m). Threat levels are:
        * No threat (H <= 0.2)
        * Marine Warning (0.2 <= H <= 0.55)
        * Land warning (H >= 0.55)
        * (Not currently supported by JATWC) Minor land warning (0.55 <= H <= 1.5)
        * (Not currently supported by JATWC) Major land warning (1.5 >= H)
    * The above thresholds are suggested in the paper by Greenslade et al. (2020) - (see https://doi.org/10.1007/s00024-019-02377-z)
    * The zones are limited to be inside the PTHA 1/10000 exceedance-rate (84th percentile epistemic uncertainty). This ensures the upper limit is not an unstable random quantity, and gives a uniform degree of conservatism at different sites. Note later on these maps were manually edited to produce the final zones (see description in AJEM manuscript), and at that point we chose to use a 1/2500 (84th pc) limit, which is more standard.

* ./elevation_contours
    * Elevation contours (0, 1, 2, ... 10 m) for ATWS zones, derived from the inundation model elevation grid. 
    * Useful reference to compare with other inundation zones.

# How to run the calculations

## Step 1
The script [compute_scenario_statistics_in_zone.R](compute_scenario_statistics_in_zone.R) is used to compute the JATWC H parameter for each scenario in a specified warning zone. For instance the computations in the "Perth Coast" zone are executed with:
```r
Rscript compute_scenario_statistics_in_zone.R 'Perth Coast'
```

## Step 2
The script [map_threat_levels_in_zone.R](map_threat_levels_in_zone.R) can make (temporary) rasters for each threat category. These rasters give the maximum waterlevel over all scenarios in that threat level category, without any consideriation of exceedance-rates from the PTHA. They represent an intermediate step in the calculations. Run with, e.g.:
```r
Rscript map_threat_levels_in_zone.R 'Perth Coast'
```

## Step 3
Finally we limit the onshore zone to a rare PTHA exceedance-rate, and then export to shapefile. The script [convert_raster_zones_to_polygons_with_PTHA_limits.R](convert_raster_zones_to_polygons_with_PTHA_limits.R) does this, with exceedance rate of 1/10000 at the 84th percentile of the epistemic uncertainty (NB: during later edits this was changed to 1/2500 @ 84th pc). It can be run for all warning types with:
```
#!/bin/bash
for warning_type in 'no_threat' 'marine_warning' 'land_warning' 'minor_land_warning' \
    'major_land_warning'; do
    Rscript convert_raster_zones_to_polygons_with_PTHA_limits.R 'Perth Coast' $warning_type;
    done;
```

## Why limit the onshore zone to a rare PTHA exceedance-rate?

* Suppose that no such limits were imposed, so the land-threat inundation zone included all areas inundated by at least one scenario. This would have two shortcomings:
    * The result would not be stable, i.e. the "all scenarios" zone might change significantly if we used a different (yet equally valid) set of random scenarios. This is because our inundation scenarios are randomly sampled from the PTHA18.
    * The "all scenarios" zone will be very conservative, but in a way that will vary from site to site, and is difficult to quantify. 

* To make the results repeatable, and explicitly control the conservatism, we can instead limit the inundation zone to a rare exceedance-rate as defined by the PTHA results. 
    * Here we used the 1/10,000 exceedance-rate for the 84th percentile epistemic uncertainty. Alternative choices could be made.
    * Visual inspection suggested that, in the Perth Coast zone, similar results would be obtained using the 1/10,000 exceedance-rate for the logic-tree-mean curve, rather than the 84th percentile, although the latter zone is often slightly larger. 
    * Our tests to date (convergence, comparison with the offshore PTHA) suggest that this exceedance-rate is reasonably well approximated with our random scenarios (i.e. we would expect similar results if we re-ran the calculations with new random scenarios).


# How to make the elevation contours
The script [make_elevation_contours.R](make_elevation_contours.R) can make elevation contours from our inundation model, for domains inside each ATWS zone. Run with (e.g.):
```r
Rscript make_elevation_contours.R "Perth Coast"
```

# How to compute the depth maxima over all marine warning scenarios

When we combine all the marine warning scenarios and take the pointwise max-stage, in some places we see minor inundation. Thus, according to the model, there could be minor inundation during some marine warning scenarios. 

To help understand how significant this inundation could be, it is useful to look at the maximum depth over all the marine warning scenarios. This can be computed with:
```r
Rscript compute_max_depths_for_marine_warning_scenarios.R
```
(Currently hard-coded for the Perth Coast zone only).

# Notes on the code

When writing this code I was trying to transition away from legacy R spatial packages (`raster`, `sp`) towards `sf` and `terra` and `stars`. But at the time of writing `terra` would not compile on NCI, and I was inexperienced with `sf` and `stars`. So various compromises were made. 

Turns out it would be possible to replace much of the `raster` and `terra` functionality with `stars`. I just didn't realise how to do it when writing. 
