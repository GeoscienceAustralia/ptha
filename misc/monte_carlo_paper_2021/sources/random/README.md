# Sample random scenarios from PTHA18
-------------------------------------

Here we sample random scenarios from PTHA18 for the Tongatapu study, and generate the scenario initial conditions. Only the `kermadectonga2` source-zone is considered. We use 'stratified+importance' sampling, where the scenarios are stratified by magnitude, and sampled based on the wave height at a site offshore of Tongatapu. This is combined with non-uniform sampling in different magnitude-bins, as specified in the file [Non_uniform_sampling_effort_compromise_stratifiedImportance.csv](Non_uniform_sampling_effort_compromise_stratifiedImportance.csv) that was created in the previous [optimal sampling calculations](../../optimal_sampling/).

The PTHA18 treats the kermadectonga2 source with a mixture of `unsegmented` and `segmented` treatments (union of 3 segments). For every case, we use THE SAME SINGLE MONTE-CARLO SAMPLE (which is sampled using weights based on the 'logic-tree-mean' rates, and the max-stage near Tongatapu). The way the code is written, it is convenient to do the calculation separately for unsegmented and each segment, but we have checks to confirm that the same sample is drawn in every case (which is achieved by re-setting the random-seed).

Importance sampling is used to put greater emphasis on scenarios with larger max-stage near Tongatapu.

The sample of scenarios is checked by plotting the max-stage exceedance-rates implied by the randomly sampled scenarios and their individual rates, at various PTHA18 gauge points. These are compared with the original PTHA18 max-stage vs exceedance-rate results. If we take a large enough sample with a good-choice of importance-sampling strategy, then these should be very similar in the range of exceedance-rates of interest.

## Key codes:

* [./select_random_scenarios.R](./select_random_scenarios.R) - This samples the scenarios (with importance sampling based on the max-stage at a site near Tongatapu), and assigns them nominal rates to achieve consistency with PTHA18. The scenarios and their rates are written to csv files (included here for convenience). We also make plots comparing the max-stage vs exceedance-rates implied by the sample of scenarios, and the original PTHA18.

* [./generate_initial_conditions.R](./generate_initial_conditions.R) - This can be run after running [select_random_scenarios.R](select_random_scenarios.R), in order to make rasters with the vertical deformation associated with each randomly selected scenario.


## How to run these codes

To run this code on NCI, start up an interactive job with access to the right filesystems (see `get_interactive_job.sh`) and do:

```r
source R_400_NCI_modules.sh

# Randomly sample scenarios, using importance sampling to preferentially pick
# those with higher PTHA18 max-stage near Tonga.
Rscript select_random_scenarios.R

# Generate initial condition tifs for the above scenarios.
# The following assumes you are running 48 cores.
Rscript generate_initial_conditions.R
```
