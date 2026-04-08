# Probabilistic hazard calculations 

## Logic tree mean inundation rate

1. Modify the logic-tree-mean exceedance rate test code to suit your case, then run it and check that all cases PASS
```
Rscript test_exceedance_rate_raster_calculations.R # Should print a few "PASS"
```
2. Modify the main logic-tree-mean exceedance rate run script to suit your case, then run it
```
qsub run_compute_exceedance_rates_for_threshold_depth_logic_tree_mean.sh
```
3. At this point I made the improved directory structure (although not everything has been created, so it will report missing files, no problem)
```    
source make_directory_structure.sh
```
4. Sum the logic-tree-mean exceedance rates over source zones
```
Rscript compute_mean_exrate_upper_CI.R ptha18-kalbarri2onslow-hazard/highres_depth_with_variance depth 0.001
```
5. Go into the new folder in `ptha18-kalbarri2onslow-hazard/highres_depth_with_variance` that contains the `sum_of_source_zones` results. Then manually make VRT files for the important ones using a command like the one below.
```
gdalbuildvrt -resolution highest __NAME_OF_VRT__.vrt string_matching_relevant_files*.tif
```

## 84th and 16th percentile exceedance rate

1. Compute exceedance rates given a variable, threshold and epistemic uncertainty. We use a script to spread the calculation over a set of qsub scripts.

```
# First check the tests work -- modify the script inputs as needed.
Rscript test_compute_exceedance_rates_at_epistemic_uncertainty.R

Rscript make_exceedance_rate_jobs.R # Modify this and the template scripts as needed

# Then submit all the resulting scripts
mkdir submitted_qsub_jobs
for i in run_compute_exceedance_rates_at_epistemic_uncertainty_depth_*.sh; do echo $i ; qsub $i; mv $i submitted_qsub_jobs/ ; done

# When they're finished, move the folders to the right location (ignore warnings about missing directories)
source make_directory_structure.sh

# And compute the sum of percentiles
Rscript compute_sum_of_percentiles.R ptha18-kalbarri2onslow-hazard/highres_depth_epistemic_uncertainty/ 84 depth 0.001
Rscript compute_sum_of_percentiles.R ptha18-kalbarri2onslow-hazard/highres_depth_epistemic_uncertainty/ 16 depth 0.001

# Then go to the folders inside
# `ptha18-kalbarri2onslow-hazard/highres_depth_epistemic_uncertainty/84pc`
# and `16pc` that contain the `sum_of_source_zones` results.
# Make VRT's in these folders using `gdalbuildvrt`
```

## Thresholds matching a given exceedance rate and epistemic uncertainty percentile

1. Compute thresholds at prescribed exceedance rates and epistemic uncertainty levels
```
# First check the test is working - modify inputs as needed
Rscript test_compute_threshold_at_exceedance_rate_of_epistemic_uncertainty.R

# If the tests PASS, then make qsub scripts to do the calculations
Rscript make_threshold_epistemic_uncertainty_jobs.R

# qsub all the scripts

# When the jobs are finished, go into the created folders and make vrt files.
```


