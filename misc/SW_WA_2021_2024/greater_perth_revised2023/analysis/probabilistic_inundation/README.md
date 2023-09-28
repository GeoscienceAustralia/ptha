# Probabilistic Inundation calculations

This folder contains scripts for probabilistic inundation calculation, which can be run after all the models in [../../swals](../../swals) have been simulated.

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

# Get the logic-tree-mean inundation exceedance-rate and variance, summed over sources.
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

# Get the inundation exceedance-rates at different epistemic uncertainty percentiles, summed over source-zones.
# This assumes co-monotonic epistemic uncertainties between the sources (conservative).
# (Runs in parallel)
Rscript compute_sum_of_percentiles.R ptha18-GreaterPerth2023-sealevel60cm/highres_depth_epistemic_uncertainty/ 84 depth 0.001
Rscript compute_sum_of_percentiles.R ptha18-GreaterPerth2023-sealevel60cm/highres_depth_epistemic_uncertainty/ 16 depth 0.001
# This created folders containing sums over source zones, with names like:
# ptha18-GreaterPerth2023-sealevel60cm/highres_depth_epistemic_uncertainty/84pc/ptha18-GreaterPerth2023-sealevel60cm-depth_exrate_0.001_0.84_sum_of_source_zones/

# At this point you should compress the output folders, as they contain many many files.

```

