# Analysis

## Check the exceedance-rate curves vs PTHA18
Start by running the tests in [probabilistic_inundation/test](probabilistic_inundation/test) and subdirectories. If these aren't good, consider doing more runs. Multiple batches can be combined using multiple importance sampling, see the [../../nsw_2023_2024](../../nsw_2023_2024) for an example of this.

## Compute the probabilistic inundation
These are split into:
- [probabilistic_inundation/exrate_given_threshold](probabilistic_inundation/exrate_given_threshold) to compute how often a threshold is exceeded. E.G. depth>0.001m gives the rate of inundation, also depending on the epistemic uncertainty treatment (logic tree mean or percentiles)
- [probabilistic_inundation/threshold_given_exrate](probabilistic_inundation/threshold_given_exrate) more arduously computes what threshold corresponds to a given exceedance rate. E.G. what's the depth for a 1/2500 exceedance rate and epistemic uncertainty percentile.

## JATWC to inundation
Using JATWC warning rules, compute the footprint covering all modelled scenarios for each current JATWC warning. 

## Arrival times
Compute the mean and minimum tsunami arrival times.

## Sea level rise
Plot the difference in results for a single scenario under a sea level rise of 0.8m. The plots compare the max stage and max speed with and without sea level rise.
