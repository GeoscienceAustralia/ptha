# Multiple importance sampling probabilistic inundation codes

The code here was modified from the single-importance-sample probabilistic inundation codes to do multiple importance sampling.

## Logic-tree mean exceedance rate and 95% confidence interval

For the mean exceedance-rate and its uncertainty, we directly compute the weighted sum (with square-weights for the variance) from previous logic-tree-mean results.

See `compute_mean_rate_from_previous_calculations.R`.

This is nothing like the single-sample code -- although it uses results of the single-sample code as input.


## 84th percentile exceedance-rate

Before running serious jobs, ensure the test passes (`test_compute_exceedance_rates_at_epistemic_uncertainty.R`).

The main code is `compute_exceedance_rates_at_epistemic_uncertainty_percentile.R`
* Note: This could do with a refactor to share more code with `compute_threshold_at_exceedance_rate_of_epistemic_uncertainty_percentile.R`. The latter necessarily puts more of the functionality into functions.

We split the work into a number of single-node jobs with `make_exceedance_rate_jobs.R`, which use the template PBS script `run_compute_exceedance_rates_at_epistemic_uncertainty_VARIABLE_SOURCEZONE_PERCENTILE_LOWER_UPPER_EXCEEDANCETHRESHOLD.sh` to create many jobs.

The associated code is quite similar to the single-sample code. But it is modified to compute the weights and sample-size associated with each sample, and integrated them into the rate calculation.

## Threshold at 1/2500 84th percentile

Before running serious jobs, ensure the test passes (`test_compute_threshold_at_exceedance_rate_of_epistemic_uncertainty.R`).

The main code is `compute_threshold_at_exceedance_rate_of_epistemic_uncertainty_percentile.R`. 

We split the work into a number of single-node jobs with `make_threshold_epistemic_uncertainty_jobs.R`. 
* This uses the template scripts `run_compute_thresholds_at_exceedance_rate_of_epistemic_uncertainty_percentile_[D-M]*.sh` that should be modified before creating the scripts to match what you need.

After running this, use the following to create new 84th percentile products where we remove values below the `uniroot` tolerance used to compute the above rasters. This is important to avoid having non-zero cells (within the uniroot tolerance) at sites that are dry at the 1/2500 84th percentile. The script will need editing first.
```
Rscript tidy_lower_bounds_in_thresholds_at_exceedance_rate_of_epistemic_uncertainty_percentile.R
``` 

Then compute "depth above initial condition at sites with elevation > 0" with something like this, after editing the script to ensure the parameter values are OK.
```
Rscript make_depth_above_initial_condition.R 
```

The associated code is quite similar to the single-sample code. But it is modified to compute the weights and sample-size associated with each sample, and integrated them into the rate calculation. It also integrates some unrelated updates which mean we don't have to pre-specify a search range for each variable.

## Threshold at 1/250 50th percentile

This is similar to the 1/2500 84th percentile case above, except the template scripts have names like `run_compute_thresholds_at_exceedance_rate_of_epistemic_uncertainty_percentile_median_1in250_[D-M]*.sh`.
* To create the run scripts, you need to modify `make_threshold_epistemic_uncertainty_jobs.R` to use these template scripts. The latter should be edited before running.

After running this, use `tidy_lower_bounds_in_thresholds_at_exceedance_rate_of_epistemic_uncertainty_percentile.R` to create new 50th percentile products where we remove values below the `uniroot` tolerance used to compute the above rasters. This is important to avoid having non-zero cells (within the uniroot tolerance) at sites that are dry at the 1/250 50th percentile.

## Flood hazard categories

When the previous calculations are done, we can create flood hazard categories with something like this (after editing the input parameters to ensure they are OK):
```
Rscript compute_flood_hazard_categories.R 
```

Notice this uses the depth above initial condition which has been limited to sites with elevation above 0. This places similar restrictions on the flood hazard category areas.

## Marine hazard thresholds

This might be best done using the threshold approach based on maximum speeds, as per Lynett et al. (2014). The thresholds are 1.5, 3, 4.5 m/s (or 3, 6, 9 knots).
