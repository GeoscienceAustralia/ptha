# Sample random scenarios from PTHA18
-------------------------------------

Here we sample random scenarios from PTHA18 for the Tongatapu study, and generate the scenario initial conditions. Only the `kermadectonga2` source-zone is considered.

The PTHA18 treats the kermadectonga2 source with a mixture of `unsegmented` and `segmented` treatments (union of 3 segments). To enable that to be carried through to the inundation PTHA, here we separately sample random scenarios for the unsegmented curves, and for each of the 3 segments. To do this we need to work with the `compute_rates_all_sources.RData` R session file, which is downloaded if we don't have it.

Importance sampling is used to put greater emphasis on scenarios with larger max-stage near Tongatapu. Each scenario is assigned an individual rate that is nominal, but leads to overall consistency with PTHA18. 

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
