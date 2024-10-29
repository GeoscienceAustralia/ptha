# Multiple importance sampling probabilistic inundation codes

The code here was modified from the single-importance-sample probabilistic inundation codes to do multiple importance sampling.

## Logic-tree mean exceedance rate and 95% confidence interval

For the mean exceedance-rate and its uncertainty, we directly compute the weighted sum (and square-weighted variance) from previous logic-tree-mean results.

See `compute_mean_rate_from_previous_calculations.R`.


## 84th percentile exceedance-rate

Before running serious jobs, ensure the test passes (`test_compute_exceedance_rates_at_epistemic_uncertainty.R`).

The main code is `compute_exceedance_rates_at_epistemic_uncertainty_percentile.R`
* Note: This could do with a refactor to share more code with `compute_threshold_at_exceedance_rate_of_epistemic_uncertainty_percentile.R`. The latter necessarily puts more of the functionality into functions.

We split the work into a number of single-node jobs with `make_exceedance_rate_jobs.R`, which use the template PBS script `run_compute_exceedance_rates_at_epistemic_uncertainty_VARIABLE_SOURCEZONE_PERCENTILE_LOWER_UPPER_EXCEEDANCETHRESHOLD.sh` to create many jobs.


## Threshold at 1/2500 84th percentile

Before running serious jobs, ensure the test passes (`test_compute_threshold_at_exceedance_rate_of_epistemic_uncertainty.R`).

The main code is `compute_threshold_at_exceedance_rate_of_epistemic_uncertainty_percentile.R`. 

We split the work into a number of single-node jobs with `make_threshold_epistemic_uncertainty_jobs.R`. 
* This uses the template scripts `run_compute_thresholds_at_exceedance_rate_of_epistemic_uncertainty_percentile_[D-M]*.sh` that should be modified before creating the scripts to match what you need.

## Threshold at 1/250 50th percentile

This is similar to the 1/2500 84th percentile case above, except the template scripts have names like `run_compute_thresholds_at_exceedance_rate_of_epistemic_uncertainty_percentile_median_1in250_[D-M]*.sh`.
* To create the run scripts, you need to modify `make_threshold_epistemic_uncertainty_jobs.R` to use these template scripts. The latter should be edited before running.
