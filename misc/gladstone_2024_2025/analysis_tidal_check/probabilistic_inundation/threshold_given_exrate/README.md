# Calculation of thresholds corresponding to an exceedance-rate and epistemic uncertainty percentile

Here we provide code to compute the depth or max-speed or max-stage or max-flux at a given exceedance-rate and epistemic uncertainty percentile (second set of code, below).
* This uses root-finding with a prescribed tolerance.

The (approximate) max-stage with a given exceedance-rate. 
  * It is approximate because we only use the max-stage values computed previously (1.1, 2.1, 3.1, ..., 10.1).
  * The computed solution is 'rounded down' from the exact solution to the nearest binned value
  * This is a simple approach to computing a quantity of interest at a given exceedance-rate.
  * For a more exact approach see the code below

```bash

# Modify compute_threshold_at_exceedance_rate_of_epistemic_uncertainty_percentile.R for your case.
# Note sub-sampling can be used to speed up the calculations (e.g. only computing the middle pixel
# in the 3x3 grid, and defining the other cells from this).

# I wrote a (non-generic) script to run highres domains only.
# This doesn't scale well so it is better to use the other scripts below
# qsub run_compute_thresholds_at_exceedance_rate_of_epistemic_uncertainty_percentile.sh
```

In general it is better to split the work into multiple jobs using the script
`make_threshold_epistemic_uncertainty_jobs.R`. 
It produces a number of separate jobs for either depth, max-stage, max-flux, max-speed 
thresholds with epistemic uncertainty. Each job will work on a subset of domains. The approach is:
```bash
# Make the submission scripts
Rscript make_threshold_epistemic_uncertainty_jobs.R
# Manually check them before proceeding.

# Submit once you're confident it's OK
for i in run_compute_thresholds_at_exceedance_rate_of_epistemic_uncertainty_percentile_[D-M]*[1-9]*.sh ;
    do echo $i ;
    qsub $i ;
    mv $i submitted_epistemic_uncertainty_jobs ;
    done
```

  * For dry regions that are wet in at least 1 scenario, [compute_threshold_at_exceedance_rate_of_epistemic_uncertainty_percentile.R](compute_threshold_at_exceedance_rate_of_epistemic_uncertainty_percentile.R) will by default return the lower bound of its search space. For max-stage, this is typically equal to the ambient sea level. In practice we want to remove this for end-user products.