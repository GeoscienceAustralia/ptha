# Probabilistic Inundation calculations

This folder contains scripts for probabilistic inundation calculation, which can be run after all the models in [../../swals](../../swals) have been simulated.

Below we show how to compute rasters depicting:
* The logic-tree-mean rate of inundation (depth > 1mm).
* The rate of inundation (depth > 1mm) at the 16th and 84th percentile epistemic uncertainty
* The logic-tree-mean exceedance-rate for a range of max-stage values (0.6, 1.6, 2.6, ...)
* The (approximate) max-stage with a given exceedance-rate. 
  * It is approximate because we only use the max-stage values computed previously (0.6, 1.6, 2.6, ...).
    * The computed solution is 'rounded down' from the exact solution to the nearest binned value
  * This is a simple approach to computing a quantity of interest at a given exceedance-rate.
    * For a more exact approach see [here](https://github.com/GeoscienceAustralia/ptha/tree/master/misc/monte_carlo_paper_2021/analysis/probabilistic_inundation).

Before running anything you'll need to modify [application_specific_inputs.R](application_specific_inputs.R) for your case.

## How to run

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
for i in run_compute_exceedance_rates_at_epistemic_uncertainty_depth*.sh; do echo $i; qsub $i; mv $i submitted_epistemic_uncertainty_jobs; done

# At this point the outputs are inside new folders within the current directory.
# It's a bit messy and should be better organised.
# Move them into organised sub-folders by editing "make_directory_structure.sh"
# for your case. Then run it.
source make_directory_structure.sh
# You should end up with a folder structure like this:
#  ptha18-GreaterPerth2023-sealevel60cm/
#
#    highres_depth_with_variance/  ## DEPTH, LOGIC TREE MEAN RESULTS
#      ptha18-GreaterPerth2023-sealevel60cm-depth-LogicTreeMean-outerrisesunda/
#      ptha18-GreaterPerth2023-sealevel60cm-depth-LogicTreeMean-sunda2/
#      ... other source zones if present ...
#
#    highres_max_stage_with_variance/  ## MAX_STAGE, LOGIC TREE MEAN RESULTS
#      ptha18-GreaterPerth2023-sealevel60cm-max_stage-LogicTreeMean-outerrisesunda/
#      ptha18-GreaterPerth2023-sealevel60cm-max_stage-LogicTreeMean-sunda2/
#      ... other source zones if present ...
#
#    highres_depth_epistemic_uncertainty/ ## DEPTH, EPISTEMIC UNCERTAINTY RESULTS
#      84pc/
#        ptha18-GreaterPerth2023-sealevel60cm-depth_exrate_0.001_0.84_outerrisesunda/
#        ptha18-GreaterPerth2023-sealevel60cm-depth_exrate_0.001_0.84_sunda2/
#        ... other source zones if present ...
#      16pc/
#        ptha18-GreaterPerth2023-sealevel60cm-depth_exrate_0.001_0.16_outerrisesunda/
#        ptha18-GreaterPerth2023-sealevel60cm-depth_exrate_0.001_0.16_sunda2/
#        ... other source zones if present ...
#
#    EMPTY FOLDERS FOR MAX STAGE EPISTEMIC UNCERTAINTIES
#

# Use calculations above to make the logic-tree-mean inundation exceedance-rate and variance, summed over sources.
# The command line argument gives the path to the logic-tree-mean results above,
# the variable (depth) and the exceedance-threshold (0.001). This uses parallel computing.
Rscript compute_mean_exrate_upper_CI.R ptha18-GreaterPerth2023-sealevel60cm/highres_depth_with_variance depth 0.001
# This created a folder inside the 'highres_depth_with_variance' sub-folder above, containing
# results summed over source-zones (exceedance-rate, upper 95% CI for true exeedance-rate, variance). 
# The folder name is like:
# ./ptha18-GreaterPerth2023-sealevel60cm/highres_depth_with_variance/ptha18-GreaterPerth2023-sealevel60cm-depth-LogicTreeMean-sum_of_source_zones

# Get the logic-tree mean exceedance-rates for all the alternative max-stage
# thresholds (as specified in the earlier qsub script). This uses parallel
# computing.
for stage_threshold in 0.601 1.6 2.6 3.6 4.6 5.6 6.6 7.6 8.6 9.6 10.6 ; do
    echo $stage_threshold ;
    Rscript compute_mean_exrate_upper_CI.R ptha18-GreaterPerth2023-sealevel60cm/highres_max_stage_with_variance max_stage $stage_threshold ;
  done

# Convert the files just created into a raster showing the max-stage threshold just
# below a given exceedance-rate (among the stage_threshold values considered
# above). This is a cheap way to make plots of the wave size matching a given
# exceedance-rate (rounded down to one of the thresholds above). 
# Example here for max_stage at an exceedance-rate of 1/500 = 0.002
Rscript compute_binned_thresholds_matching_exrate_from_set_of_exrate_rasters.R ptha18-GreaterPerth2023-sealevel60cm/highres_max_stage_with_variance/ptha18-GreaterPerth2023-sealevel60cm-max_stage-LogicTreeMean-sum_of_source_zones max_stage 0.002

# Get the inundation exceedance-rates at different epistemic uncertainty percentiles, summed over source-zones.
# This assumes co-monotonic epistemic uncertainties between the sources (conservative).
# (Runs in parallel)
Rscript compute_sum_of_percentiles.R ptha18-GreaterPerth2023-sealevel60cm/highres_depth_epistemic_uncertainty/ 84 depth 0.001
Rscript compute_sum_of_percentiles.R ptha18-GreaterPerth2023-sealevel60cm/highres_depth_epistemic_uncertainty/ 16 depth 0.001
# This created folders containing sums over source zones, with names like:
# ptha18-GreaterPerth2023-sealevel60cm/highres_depth_epistemic_uncertainty/84pc/ptha18-GreaterPerth2023-sealevel60cm-depth_exrate_0.001_0.84_sum_of_source_zones/

# At this point (or before) you should compress the output folders, as they can
# contain many many files (e.g. 10^5 for Greater Perth). 

```

