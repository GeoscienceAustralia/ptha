# Sample scenarios from a set of source zones using basic importance sampling. 

The total number of scenarios on each source zone was defined subjectively, by examining PTHA18 deaggregations and including sources that might lead to large tsunamis in eastern Australia. Smaller sources are given stronger representation than might be suggested by PTHA18, because it is possible that they show different onshore behaviour (not covered by PTHA18).

A difference with previous work is that we don't use stratified sampling by magnitude bins. On some sources we would not have enough scenarios to sample each PTHA18 magnitude. Thus we would have to merge bins or dispense with them entirely. 

Herein we ignore magnitude bins. Instead we sample scenarios to ensure roughly equal counts between the following max-stage values at the target point: 
* `0.0, STG250_84, STG2500_84, STG_10000_84, MAX_STAGE_FROM_SCENARIO`. 
  * The notation `STG250_84` means the max-stage at the 84th percentile 1/250 exceedance-rate at the target point.
* Additional stage-bin boundaries are defined by interpolating between those above.
  * This further promotes even sampling of stages 
* The method is implemented using importance sampling. 
  * The number of scenarios in each bin is random, not fixed as in stratified sampling.
* An alternative that was considered but not implemented was to use evenly spaced stage bin boundaries 
  * Such as `0, 0.1, 0.2, 0.3, ... MAX_STAGE_ROUNDED_UP`.
  * This wasn't done because in our application, the return periods of interest will be more concentrated around particular return periods. 
    * Even bin boundaries could easily lead to heavy sampling of very rare scenarios.

## How to run it
First create a directory for you new batch. Copy a `sampling_config.R` file in and adjust as required.
Use the `run_create_scenarios.pbs` script which does the following:
```bash
source R_431_NCI_modules.sh # On NCI only
cd batch_1  # change to your directory name of your batch
Rscript ../create_scenarios.R
cd ..
qsub run_create_scenarios.pbs  # edit
```

## Key codes

* sampling_config.R - Edit this to control the sampling approach
* importance_sampling_utilities.R - Useful functions for importance sampling
* create_scenarios.R - Main script for sampling scenarios and making test plots.
* create_initial_conditions_for_scenarios.R - Make initial condition rasters once scenarios have been created.
