# Probabilistic hazard calculations 

## Logic tree mean inundation rate

1. Modify the logic-tree-mean exceedance rate test code to suit your case, then run it and check that all cases PASS
```{bash}
Rscript test_exceedance_rate_raster_calculations.R # Should print a few "PASS"
```
2. Modify the main logic-tree-mean exceedance rate run script to suit your case, then run it
```{bash}
qsub run_compute_exceedance_rates_for_threshold_depth_logic_tree_mean.sh
```
3. At this point I made the improved directory structure (although not everything has been created, so it will report missing files, no problem)
```{bash} 
source make_directory_structure.sh
```
4. Sum the logic-tree-mean exceedance rates over source zones (uses parallel computing)
```{bash}
Rscript compute_mean_exrate_upper_CI.R ptha18-coralbay2onslow-hazard/highres_depth_with_variance depth 0.001
```
5. Go into the new folder in `ptha18-coralbay2onslow-hazard/highres_depth_with_variance` that contains the `sum_of_source_zones` results. Then manually make VRT files for the important ones using a command like the one below:
```{bash}
gdalbuildvrt -resolution highest __NAME_OF_VRT__.vrt string_matching_relevant_files*.tif
```

## 84th and 16th percentile exceedance rate

1. Compute exceedance rates given a variable, threshold and epistemic uncertainty. We use a script to spread the calculation over a set of qsub scripts.
First check the tests work -- modify the script inputs as needed.

```{bash}
Rscript test_compute_exceedance_rates_at_epistemic_uncertainty.R

Rscript make_exceedance_rate_jobs.R # Modify this and the template scripts as needed
```

2. Then submit all the resulting scripts
```{bash}
for i in exrate_given_epi_depth_*.pbs; do echo $i ; qsub $i; mv $i log/ ; done
```
When they're finished, move the folders to the right location (ignore warnings about missing directories)
```{bash}
. make_directory_structure.sh
```

3. And compute the sum of percentiles
```{bash}
Rscript compute_sum_of_percentiles.R ptha18-coralbay2onslow-hazard/highres_depth_epistemic_uncertainty/ 84 depth 0.001
Rscript compute_sum_of_percentiles.R ptha18-coralbay2onslow-hazard/highres_depth_epistemic_uncertainty/ 16 depth 0.001
```

4. Then go to the folders inside
`ptha18-coralbay2onslow-hazard/highres_depth_epistemic_uncertainty/84pc`
and `16pc` that contain the `sum_of_source_zones` results.
Make VRT's in these folders using `gdalbuildvrt`.

## Thresholds matching a given exceedance rate and epistemic uncertainty percentile

1. Compute thresholds at prescribed exceedance rates and epistemic uncertainty levels

First check the test is working - modify inputs as needed
```{bash}
Rscript test_compute_threshold_at_exceedance_rate_of_epistemic_uncertainty.R
```

2. If the tests PASS, then make qsub scripts to do the calculations
```
Rscript make_threshold_epistemic_uncertainty_jobs.R
```

3. qsub all the scripts as above:
```{bash}
for i in threshold_given_exrate_epi_percentile_*.pbs; do echo $i ; qsub $i; mv $i log/ ; done
```

4. Some of the jobs might fail due to running out of time. 
You could make new qsub scripts for the cases that didn't complete by suitably modifying
```{bash}
    make_failed_threshold_epistemic_uncertainty_jobs.R
```
Alternatively, consider manually editing the failed qsub scripts (e.g. to increase the time)

5. When the jobs are finished, go into the created folders and make vrt files.
