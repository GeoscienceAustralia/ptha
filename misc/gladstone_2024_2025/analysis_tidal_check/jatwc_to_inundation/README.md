# Compute inundation zones corresponding to JATWC warning categories, with upper limits based on PTHA.

Once all scripts have been run, the key outputs are polygon shapefiles in folders named like `Inundation_zones/{zone_name}/{zone_name}_{warning_type}_with_PTHA-exrate-limit_{...}/`

# Details of folder contents (once scripts have been run)

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
    * The zones are limited to be inside the PTHA 1/2500 exceedance-rate (84th percentile epistemic uncertainty). This ensures the upper limit is not an unstable random quantity, and gives a uniform degree of conservatism at different sites.

* ./elevation_contours
    * Elevation contours (0, 1, 2, ... 10 m) for ATWS zones, derived from the inundation model elevation grid. 
    * Useful reference to compare with other inundation zones.

* ./domains_shapefile
    * Shapefile showing the locations of domains in the model. For outputs that consist of sets of rasters, this provides an overview of where each raster is located. It was created using the script [../../swals/make_domains_shapefile.R](../../swals/make_domains_shapefile.R) and manually copied here (for convenience only).

# How to run the calculations

## Step 0
Copy one example of the elevation raster files to a local directory using [copy_elevation_rasters_locally.R](copy_elevation_rasters_locally.R). The folder from which we copy is specified in the latter script, and must have exactly the same domain structure as all other hazard scenario models.
```r
Rscript copy_elevation_rasters_locally.R
```

Edit [application_specific_metadata.R](application_specific_metadata.R) to match the model setup.

## Step 1
The script [compute_scenario_statistics_in_zone.R](compute_scenario_statistics_in_zone.R) is used to compute the JATWC H parameter for each scenario in a specified warning zone. For instance the computations in the "Capricornia Coast" zone are executed with:
```r
Rscript compute_scenario_statistics_in_zone.R 'Capricornia Coast'
```

## Step 2
The script [map_threat_levels_in_zone.R](map_threat_levels_in_zone.R) can make (temporary) rasters for each threat category. These rasters give the maximum waterlevel over all scenarios in that threat level category, without any consideriation of exceedance-rates from the PTHA. They represent an intermediate step in the calculations. Run with, e.g.:
```r
Rscript map_threat_levels_in_zone.R 'Capricornia Coast'
```

## Step 3
Finally we limit the onshore zone to a rare PTHA exceedance-rate, and then export to shapefile. The script [convert_raster_zones_to_polygons_with_PTHA_limits.R](convert_raster_zones_to_polygons_with_PTHA_limits.R) does this, with exceedance rate of 1/2500 at the 84th percentile of the epistemic uncertainty. It can be run for all warning types with:
```
#!/bin/bash 
for warning_type in 'no_threat' 'marine_warning' 'land_warning' 'minor_land_warning' 'major_land_warning'; do
    Rscript convert_raster_zones_to_polygons_with_PTHA_limits.R 'Capricornia Coast' $warning_type;
    done;
```

### Issues with flooding of low-lying areas

The modelling uses a HAT sea-level at the open coast, which can flood some
usually-dry areas in estuaries where the tide is attenuated.

One way to recognize this is to look for sites in estuaries that are inside the no-threat polygon.

To recognize these sites, one useful product is the
`depth_above_initial_condition` created in the
[../probabilistic_inundation](probabilistic_inundation) folder. At low-lying
sites that only have a small modelled tsunami, this will be small.

## Why limit the onshore zone to a rare PTHA exceedance-rate?

* Suppose that no such limits were imposed, so the land-threat inundation zone included all areas inundated by at least one scenario. This would have two shortcomings:
    * The result would not be stable, i.e. the "all scenarios" zone might change significantly if we used a different (yet equally valid) set of random scenarios. This is because our inundation scenarios are randomly sampled from the PTHA18.
    * The "all scenarios" zone will be very conservative, but in a way that will vary from site to site, and is difficult to quantify. 

* To make the results repeatable, and explicitly control the conservatism, we can instead limit the inundation zone to a rare exceedance-rate as defined by the PTHA results. 
    * Here we used the 1/2500 exceedance-rate for the 84th percentile epistemic uncertainty. 
    * At the time of writing (first stage of NSW results) it isn't completely clear how well resolved this percentile is -- check that in future work.


# How to compute the depth maxima over all marine warning scenarios

When we combine all the marine warning scenarios and take the pointwise max-stage, in some places we see minor inundation. Thus, according to the model, there could be minor inundation during some marine warning scenarios. 

To help understand how significant this inundation could be, it is useful to look at the maximum depth over all the marine warning scenarios. This can be computed with (e.g.):
```r
Rscript compute_max_depths_for_marine_warning_scenarios.R "Capricornia Coast"
```

# Notes on the code

When writing this code I was trying to transition away from legacy R spatial packages (`raster`, `sp`) towards `sf` and `terra` and `stars`. But at the time of writing `terra` would not compile on NCI, and I was inexperienced with `sf` and `stars`. So various compromises were made. 

Turns out it would be possible to replace much of the `raster` and `terra` functionality with `stars`. I just didn't realise how to do it when writing. 

Subsequently I got both `stars` and `terra` working fine on NCI ([the key step involved some changes to the s2 package, see here](https://github.com/r-spatial/s2/issues/199)). But I haven't totally updated these codes. If you try to do that, beware that subtle differences can arise in interpreting NA values in rasters produced by `raster`, `terra` and `stars`, which in some situations can cause code to break. Double check for this if you do re-write.
