# Probabilistic Inundation calculations

This folder contains scripts for probabilistic inundation calculation, which can be run after all the models in [../../swals](../../swals) have been simulated.

Before running anything you'll need to modify [application_specific_inputs.R](application_specific_inputs.R) for your case.

You will also have to hand-modify various scripts below. This is explained in comments but perhaps not comprehensively. It is always a good idea to read any script prior to submitting it.


# Exceedance-rate calculations

All the calculations in this section are based around exceedance-rates. Below we show how to compute rasters depicting:
* The logic-tree-mean rate of inundation (depth > 1mm).
* The rate of inundation (depth > 1mm) at the 16th and 84th percentile epistemic uncertainty
* The logic-tree-mean exceedance-rate for a range of max-stage values (0.6, 1.6, 2.6, ...)
* The (approximate) max-stage with a given exceedance-rate. 
  * It is approximate because we only use the max-stage values computed previously (0.6, 1.6, 2.6, ...).
    * The computed solution is 'rounded down' from the exact solution to the nearest binned value
  * This is a simple approach to computing a quantity of interest at a given exceedance-rate.
    * For a more exact approach see the next section (which applies to percentile uncertainties)

```bash

# Modify application_specific_inputs.R for your case

# Inundation exceedance-rate calculations (logic-tree-mean case) 1st step
#
# Modify the test code to suit your case, then run it and check that all cases PASS
Rscript test_exceedance_rate_raster_calculations.R # Should print a few "PASS"
# Modify the main run script to suit your case, then run it
qsub run_compute_exceedance_rates_for_threshold_depth_logic_tree_mean.sh

# Exceedance-rate calculations for a range of max-stage values (logic-tree-mean case, 1st step)
# Modify the main run script to suit your case, then run it
qsub run_compute_exceedance_rates_for_threshold_max_stages_logic_tree_mean.sh

# Inundation exceedance-rate calculations (epistemic uncertainty case) 1st step
#
# Modify the test code to suit your case then run and check that all checks PASS 
# This uses mutliple cores so probably needs an interactive job on NCI.
Rscript test_compute_exceedance_rates_at_epistemic_uncertainty.R # Should print "PASS"
# Check and edit make_exceedance_rate_jobs.R and the associated template script.
# These will make a set of qsub files to run the calculations
# Make the qsub files with:
Rscript make_exceedance_rate_jobs.R
# Check a few of those qsub files (in case you made a mistake in setup)
# If they are OK then submit them, moving to a folder after submission.
mkdir -p submitted_epistemic_uncertainty_jobs
for i in run_compute_exceedance_rates_at_epistemic_uncertainty_depth*.pbs; do echo $i; qsub $i; mv $i submitted_epistemic_uncertainty_jobs; done
# Now alter the script `make_exceedance_rate_jobs.R` to the other variable: max_stage or depth and repeat the above for .. qsub; done loop:
for i in run_compute_exceedance_rates_at_epistemic_uncertainty_max_stage*.pbs; do echo $i; qsub $i; mv $i submitted_epistemic_uncertainty_jobs; done

# At this point the outputs are inside new folders within the current directory.
# It's a bit messy and should be better organised.
# Move them into organised sub-folders by editing "make_directory_structure.sh"
# for your case. Then run it.
source make_directory_structure.sh
# You should end up with a folder structure like this:
#  ptha18-midwest-sealevel60cm/
#
#    highres_depth_with_variance/  ## DEPTH, LOGIC TREE MEAN RESULTS
#      ptha18-midwest-sealevel60cm-depth-LogicTreeMean-outerrisesunda/
#      ptha18-midwest-sealevel60cm-depth-LogicTreeMean-sunda2/
#      ... other source zones if present ...
#
#    highres_max_stage_with_variance/  ## MAX_STAGE, LOGIC TREE MEAN RESULTS
#      ptha18-midwest-sealevel60cm-max_stage-LogicTreeMean-outerrisesunda/
#      ptha18-midwest-sealevel60cm-max_stage-LogicTreeMean-sunda2/
#      ... other source zones if present ...
#
#    highres_depth_epistemic_uncertainty/ ## DEPTH, EPISTEMIC UNCERTAINTY RESULTS
#      84pc/
#        ptha18-midwest-sealevel60cm-depth_exrate_0.001_0.84_outerrisesunda/
#        ptha18-midwest-sealevel60cm-depth_exrate_0.001_0.84_sunda2/
#        ... other source zones if present ...
#      16pc/
#        ptha18-midwest-sealevel60cm-depth_exrate_0.001_0.16_outerrisesunda/
#        ptha18-midwest-sealevel60cm-depth_exrate_0.001_0.16_sunda2/
#        ... other source zones if present ...
#
#    EMPTY FOLDERS FOR MAX STAGE EPISTEMIC UNCERTAINTIES
#

# Now get the logic-tree-mean inundation exceedance-rate and variance, summed over sources.
# The command line argument gives the path to the logic-tree-mean results above,
# the variable (depth) and the exceedance-threshold (0.001). This uses parallel computing.
Rscript compute_mean_exrate_upper_CI.R ptha18-midwest-sealevel60cm/highres_depth_with_variance depth 0.001
# This created a folder inside the 'highres_depth_with_variance' sub-folder above, containing
# results summed over source-zones (exceedance-rate, upper 95% CI for true exeedance-rate, variance). 
# The folder name is like:
# ./ptha18-midwest-sealevel60cm/highres_depth_with_variance/ptha18-midwest-sealevel60cm_of_source_zones

# Use calculations above to get the logic-tree mean exceedance-rates for all
# the alternative max-stage thresholds. The thresholds must be the same as
# specified in the earlier qsub script. This uses parallel computing.
for stage_threshold in 0.601 1.6 2.6 3.6 4.6 5.6 6.6 7.6 8.6 9.6 10.6 ; do
    echo $stage_threshold ;
    Rscript compute_mean_exrate_upper_CI.R ptha18-midwest-sealevel60cm/highres_max_stage_with_variance max_stage $stage_threshold ;
done

# Convert the files just created into a raster showing the max-stage threshold just
# below a given exceedance-rate (among the stage_threshold values considered
# above). This is a cheap way to make plots of the wave size matching a given
# exceedance-rate (rounded down to one of the thresholds above). 
# Example here for max_stage at an exceedance-rate of 1/500 = 0.002
Rscript compute_binned_thresholds_matching_exrate_from_set_of_exrate_rasters.R ptha18-midwest-sealevel60cm/highres_max_stage_with_variance/ptha18-midwest-sealevel60cm-max_stage-LogicTreeMean-sum_of_source_zones max_stage 0.002

# Use calculations above to get the inundation exceedance-rates at different
# epistemic uncertainty percentiles, summed over source-zones. This assumes
# co-monotonic epistemic uncertainties between the sources (conservative).
# (Runs in parallel)
Rscript compute_sum_of_percentiles.R ptha18-midwest-sealevel60cm/highres_depth_epistemic_uncertainty/ 84 depth 0.001
Rscript compute_sum_of_percentiles.R ptha18-midwest-sealevel60cm/highres_depth_epistemic_uncertainty/ 16 depth 0.001
# This created folders containing sums over source zones, with names like:
# ptha18-midwest-sealevel60cm/highres_depth_epistemic_uncertainty/84pc/ptha18-midwest-sealevel60cm-depth_exrate_0.001_0.84_sum_of_source_zones/

# At this point (or before) you should compress the output folders, as they can
# contain many many files (e.g. 10^5 for Greater Perth). 

```

# Calculation of depth / max-stage / max-flux / max-speed at a given exceedance-rate and epistemic uncertainty percentile.

Here we show how to calculate the depth / max-stage / max-flux / max-speed at a given exceedance-rate and epistemic uncertainty percentile.
* Root-finding is used to compute the variable of interest within a prescribed tolerance.

```bash

# Modify compute_threshold_at_exceedance_rate_of_epistemic_uncertainty_percentile.R for your case.
# Note sub-sampling can be used to speed up the calculations (e.g. only computing the middle pixel
# in the 3x3 grid, and defining the other cells from this).

# Modify the test code below to suit your case, then run it and check that it prints PASS.
# This will require a node, with as many cores needed in application_specific_inputs.R::MC_CORES_SR.
Rscript test_compute_threshold_at_exceedance_rate_of_epistemic_uncertainty.R

# If the test works, proceed with calculations of interest.
# I wrote a (non-generic) script to run highres domains only
#   qsub run_compute_thresholds_at_exceedance_rate_of_epistemic_uncertainty_percentile.sh
# However this tends to have problems with compute time and memory
#
# BELOW WE SHORE A MORE GENERIC AND PRACTICAL APPROACH
# We compute depth/max-stage/max-flux/max-speed at a given percentile and exceedance-rate.
# Because the calculations can be heavy we distribute them among multiple jobs. 
#
# First have a look at the template scripts named 
#  run_compute_thresholds_at_exceedance_rate_of_epistemic_uncertainty_percentile_DEPTH___MIN_DOMAIN_INDEX_____MAX_DOMAIN_INDEX__.sh
#  run_compute_thresholds_at_exceedance_rate_of_epistemic_uncertainty_percentile_MAX_FLUX___MIN_DOMAIN_INDEX_____MAX_DOMAIN_INDEX__.sh
#  run_compute_thresholds_at_exceedance_rate_of_epistemic_uncertainty_percentile_MAX_SPEED___MIN_DOMAIN_INDEX_____MAX_DOMAIN_INDEX__.sh
#  run_compute_thresholds_at_exceedance_rate_of_epistemic_uncertainty_percentile_MAX_STAGE___MIN_DOMAIN_INDEX_____MAX_DOMAIN_INDEX__.sh
# and the script
#   make_threshold_epistemic_uncertainty_jobs.R
# 
# The latter script is used to generate many scripts using the templates. Make
# sure you are happy with the template scripts and the variables in
# make_threshold_epistemic_uncertainty_jobs.R and then do
#
Rscript make_threshold_epistemic_uncertainty_jobs.R
#
# to make a bunch of PBS job scripts that run everything. This will not submit
# any calculations. But the resulting set of PBS job scripts will each do
# calculations for subset of domains for one particular flow variable. The
# number of jobs scripts per flow variable is specified inside the R script,
# and might be edited depending on how much work there is.
#
# Before submitting those scripts, eyeball them to ensure they are doing what you want.
# Then submit them
# 

# A shortcoming of the calculations above is that dry areas are not
# consistently set to NA. For depth and max-stage, areas that are flooded by at
# least one scenario but are dry at the specified exceedance-rate/percentile
# will be set to the lower-limit of the root-finding space. For max-speed and
# max-flux, even "always dry" areas will get zero values. This is easy to
# misinterpret, so we use a separate script to clean up the wet-dry values.
# 
# Make sure you edit it to set the right folders/thresholds etc, then do
Rscript tidy_lower_bounds_in_thresholds_at_exceedance_rate_of_epistemic_uncertainty_percentile.R
# I suggest using thresholds slightly larger than the lower-bound plus the tolerance.
```

# Calculation of arrival time minima and average for each source-zone

* The arrival time (minimum and average) in seconds post-earthquake, for each source-zone.
  * Our SWALS calculations define the arrival time as the minimum time at which `max_stage > (0.01 + background_stage)` and the cell is wet. The latter constraint only matters on land where the elevation exceeds the background stage.
  * The average arrival time is computed as a naive average over all modelled scenarios (ignoring scenario rates). 
    * If a scenario has cells with NA arrival times, then for those cells, the scenario is dropped from the average. This can happen because the tsunami does not exceed the threshold, or because the site is not inundated.
  * We choose to ignore scenario rates when computing the average. If the average were weighted by the scenario rates, then the results would be dominated by a few small scenarios.
  * Arrival times as defined here have the following potentially surprising properties.
    * The arrival times can have discontinuities (e.g. if a wave slightly exceeds the threshold at some sites and is slightly below at others). These propagate through to the minima and average arrival times.
    * It is possible for the average arrival time on land to be earlier than the average arrival time closer to the ocean. 
      * For example, suppose only one scenario floods the landward site and has a relatively early arrival time, while sites closer to the coast are flooded by a range of scenarios including those with much later arrival times.
      * This will not happen with minimum arrival times.

```
# Modify compute_arrival_time_minima_and_scenario_average.R 
# to fit your case.
#
# Then run with
qsub run_compute_arrival_time_minima_and_scenario_average.sh
#
```
