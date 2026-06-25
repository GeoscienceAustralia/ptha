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
    * For each modelled scenario, in each warning zone, the threat level is defined with a statistic H (representing the 95th percentile of the tsunami maximum stage in the zone, where depths > 20m). Threat levels (excepting offshore islands) are:
        * No threat (H <= 0.2)
        * Marine Warning (0.2 <= H <= 0.55)
        * Land warning (H >= 0.55)
        * (Not currently supported by JATWC) Minor land warning (0.55 <= H <= 1.5)
        * (Not currently supported by JATWC) Major land warning (1.5 >= H)
        * At offshore island ATWS zones instead use a Marine warning lower limit of 0.1, and land warning lower limit of 0.5.
            * Cocos Island, Christmas Island, Willis Island, Norfolk Island, Lord Howe Island
    * The above thresholds are suggested in the paper by Greenslade et al. (2020) - (see https://doi.org/10.1007/s00024-019-02377-z)
        * Except for the offshore island thresholds, which are in "Greenslade and Allen (2019) On the optimal amplitude thresholds for tsunami warning, Bureau Research Report - 039" and have been confirmed by discussion with JATWC personnel.
    * The zones are limited to be inside a rare PTHA exceedance-rate (84th percentile epistemic uncertainty). This ensures the upper limit is not an unstable random quantity, and gives a uniform degree of conservatism at different sites. (Initially we tried 1/10000 84th percentile but ultimately settled on 1/2500 84th percentile).

* ./elevation_contours
    * Elevation contours (0, 1, 2, ... 10 m) for ATWS zones, derived from the inundation model elevation grid WITHOUT TIDAL ADJUSTMENT (i.e. "real" elevations). 
    * Useful reference to compare with other inundation zones.

* ./elevation_in_model_no_tidal_adjustment
    * elevation as seen by the model, without a tidal adjustment

* ./elevation_in_model_with_tidal_adjustment
    * elevation as seen by the model, including a spatially varying vertical adjustment so that the initial shoreline is an fairly high tide (maximum from tidal model TPXO9v5a).

* ./domains_shapefile
    * Shapefile showing the locations of domains in the model. For outputs that consist of sets of rasters, this provides an overview of where each raster is located. It was created using the script [../../swals/make_domains_shapefile.R](../../swals/make_domains_shapefile.R) and manually copied here (for convenience only).

* ./tidal_adjustment
    * Raster showing the difference between the elevation with/without tidal adjustment. 
    * In most places this is identical to the tidal offset used. 
        * There could be exceptions in a few pixels on the coarsest domains, near wet/dry areas, at sites that we do not anyway use for inundation modelling. 
            * That's because SWALS adjusts the elevation on such domains, at sites below the ambient sea level, to be at least 1 m below the sea level AFTER tidal adjustment is applied. 
            * This is done because the coarse global domains use linear solvers with no wetting and drying, and it can help them nest with nonlinear domains without making negative depths.

# How to run the calculations

## Step 0
Copy one example of the elevation raster files to a local directory using [copy_elevation_rasters_locally.R](copy_elevation_rasters_locally.R). The folder(s) from which we copy is specified in the latter script, and must have exactly the same domain structure and elevation data as all other hazard scenario models (we get versions with/without tidal adjustment).
```r
Rscript copy_elevation_rasters_locally.R
```

Then go inside the [tidal_adjustment](tidal_adjustment) directory and do
```r
Rscript get_tidal_adjustment.R
```

Then copy information on the elevation source data here
```r
Rscript copy_elevation_source_file_index_locally.R
```

Also copy a `domains_shapefile` to the current directory. The required files were probably already created using [../../swals/make_domains_shapefile.R](../../swals/make_domains_shapefile.R). This is not necessary for the workflows below, but is convenient.

Then edit [application_specific_metadata.R](application_specific_metadata.R) to match the model setup.

## Step 1
The script [compute_scenario_statistics_in_zone.R](compute_scenario_statistics_in_zone.R) is used to compute the JATWC H parameter for each scenario in a specified warning zone. For instance: 
```r
Rscript compute_scenario_statistics_in_zone.R 'Geraldton Coast'
Rscript compute_scenario_statistics_in_zone.R 'Gascoyne Coast'
Rscript compute_scenario_statistics_in_zone.R 'Ningaloo Coast'
# (Skip Pilbara since this model only reaches Coral Bay in mid Ningaloo)
# Rscript compute_scenario_statistics_in_zone.R 'Pilbara Coast West'

```

## Step 2
The script [map_threat_levels_in_zone.R](map_threat_levels_in_zone.R) can make (temporary) rasters for each threat category. These rasters give the maximum water-level over all scenarios in that threat level category, without any consideration of exceedance-rates from the PTHA. They represent an intermediate step in the calculations. Run with, e.g.:
```r
Rscript map_threat_levels_in_zone.R 'Geraldton Coast'
Rscript map_threat_levels_in_zone.R 'Gascoyne Coast'
Rscript map_threat_levels_in_zone.R 'Ningaloo Coast'
# (Skip Pilbara since this model only reaches Coral Bay in mid Ningaloo)
# Rscript map_threat_levels_in_zone.R 'Pilbara Coast West'

```

## Step 3
Finally we limit the onshore zone to a rare PTHA exceedance-rate, and then export to shapefile. The script [convert_raster_zones_to_polygons_with_PTHA_limits.R](convert_raster_zones_to_polygons_with_PTHA_limits.R) does this. It can be run for all warning types with:
```
qsub run_all_zones.sh
```

This includes a special treatment of areas where the initial sea level was set lower than normal.

## Why limit the onshore zone to a rare PTHA exceedance-rate?

* Suppose that no such limits were imposed, so the land-threat inundation zone included all areas inundated by at least one scenario. This would have two shortcomings:
    * The result would not be stable, i.e. the "all scenarios" zone might change significantly if we used a different (yet equally valid) set of random scenarios. This is because our inundation scenarios are randomly sampled from the PTHA18.
    * The "all scenarios" zone will be very conservative, but in a way that will vary from site to site, and is difficult to quantify. 

* To make the results repeatable, and explicitly control the conservatism, we can instead limit the inundation zone to a rare exceedance-rate as defined by the PTHA results. 
    * Here we used the 1/2500 exceedance-rate for the 84th percentile epistemic uncertainty, after some initial experimentation with 1/10000 84th percentile.
    * Our tests to date (convergence, comparison with the offshore PTHA) suggest that this exceedance-rate is reasonably well approximated with our random scenarios (i.e. we would expect similar results if we re-ran the calculations with new random scenarios).


# How to make the elevation contours
The script [make_elevation_contours.R](make_elevation_contours.R) can make elevation contours from our inundation model, for domains inside each ATWS zone. Run with (e.g.):
```r
Rscript make_elevation_contours.R "Geraldton Coast"
```

# How to compute the depth maxima over all marine warning scenarios

When we combine all the marine warning scenarios and take the point-wise max-stage, in some places we see minor inundation. Thus, according to the model, there could be minor inundation during some marine warning scenarios. 

To help understand how significant this inundation could be, it is useful to look at the maximum depth over all the marine warning scenarios. This can be computed with (e.g.):
```r
Rscript compute_max_depths_for_marine_warning_scenarios.R "Geraldton Coast"
```
Beware that the rasters produced with this script apply only to the given coastal zone, and will be discontinuous with the rasters produced for neighbouring coastal zones. For example, the "Lancelin Coast" results will overlap into the "Geraldton Coast" zone, but the former are derived for scenarios having "Lancelin Coast Marine Warning", which is a different set of scenarios than those producing "Geraldton Coast Marine Warning".

# Notes on the code

When writing this code I was trying to transition away from legacy R spatial packages (`raster`, `sp`) towards `sf` and `terra` and `stars`. But at the time of writing `terra` would not compile on NCI, and I was inexperienced with `sf` and `stars`. So various compromises were made. 

Turns out it would be possible to replace much of the `raster` and `terra` functionality with `stars`. I just didn't realise how to do it when writing. 

Subsequently I got both `stars` and `terra` working fine on NCI ([the key step involved some changes to the s2 package, see here](https://github.com/r-spatial/s2/issues/199)). But I haven't totally updated these codes. If you try to do that, beware that subtle differences can arise in interpreting NA values in rasters produced by `raster`, `terra` and `stars`, which in some situations can cause code to break. Double check for this if you do re-write.
