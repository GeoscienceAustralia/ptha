# Calculation of thresholds corresponding to an exceedance-rate and epistemic uncertainty percentile

Here we provide code to compute the depth or max-speed or max-stage or max-flux at a given exceedance-rate and epistemic uncertainty percentile (second set of code, below).
* This uses root-finding with a prescribed tolerance.

## Binned thresholds (no longer recommended)
The (approximate) max-stage with a given exceedance-rate. 
  * It is approximate because we only use the max-stage values computed previously (1.1, 2.1, 3.1, ..., 10.1).
  * The computed solution is 'rounded down' from the exact solution to the nearest binned value
  * This is a simple approach to computing a quantity of interest at a given exceedance-rate.
  * For a more exact approach see the code below

## Uniroot to find hazard values at design levels
Modify compute_threshold_at_exceedance_rate_of_epistemic_uncertainty_percentile.R for your case.
Note sub-sampling can be used to speed up the calculations (e.g. only computing the middle pixel
in the 3x3 grid, and defining the other cells from this).

I wrote a (non-generic) script to run highres domains only.
This doesn't scale well so it is better to use the other scripts below
```bash
# qsub run_compute_thresholds_at_exceedance_rate_of_epistemic_uncertainty_percentile.sh
```

In general it is better to split the work into multiple jobs using the script
`make_threshold_epistemic_uncertainty_jobs.R`. 
It produces a number of separate jobs for either depth, max-stage, max-flux, max-speed 
thresholds with epistemic uncertainty. Each job will work on a subset of domains. The approach is:
1. Make the submission scripts
```bash
Rscript make_threshold_epistemic_uncertainty_jobs.R
```

2. Manually check them before proceeding. Note it takes significant compute.

3. Submit once you're confident it's OK
```bash
for pbs_file in run_compute_thresholds_at_exceedance_rate_of_epistemic_uncertainty_percentile_m*.pbs; do
  echo $pbs_file
  qsub $pbs_file
  mv $pbs_file submitted_jobs
done
```

- For dry regions that are wet in at least 1 scenario, [compute_threshold_at_exceedance_rate_of_epistemic_uncertainty_percentile.R](compute_threshold_at_exceedance_rate_of_epistemic_uncertainty_percentile.R) will by default return the lower bound of its search space. For max-stage, this is typically equal to the ambient sea level. In practice we want to remove this for end-user products.

4. Move only the successful logs, which can be useful to see if any of the jobs failed:
```bash
files=$(find . -type f -name "*.pbs.o*")
for file in $files
do
    # if it has "Exit Status:        0" mv to log/
    if grep -q "Exit Status:        0" $file
    then
        echo "Moving $file to log/"
        mv $file log/
        # and move the corresponding .pbs.e* file
        mv ${file%.o*}.e* log/
    fi
done
```

Or simply move all the logs `mv *.pbs.e* log/` and `mv *.pbs.o* log/`
