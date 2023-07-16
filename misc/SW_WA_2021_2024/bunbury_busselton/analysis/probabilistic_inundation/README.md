# Probabilistic Inundation calculations

This folder contains scripts for probabilistic inundation calculation, which can be run after all the models in [../../swals](../../swals) have been simulated.

## Key outputs (after all scripts have been run)

The main output files are in the folder [ptha18-BunburyBusseltonRevised-sealevel60cm](ptha18-BunburyBusseltonRevised-sealevel60cm). They show the exceedance-rates of inundation (logic tree mean and percentiles), as well as exceedance-rates of the max-stage above the assumed background sea-level. See the README in that folder for more information.

## Key files

* [compute_exceedance_rates_for_threshold_depth_logic_tree_mean_newParallelPartition.R](compute_exceedance_rates_for_threshold_depth_logic_tree_mean_newParallelPartition.R) - This estimates the logic-tree mean exceedance-rate for a specified depth threshold, as well as the variance of the Monte-Carlo error. It only operates on one source-zone at a time, so results of source-zones must be later summed.
    * The computational core underlying these calculations is in [exceedance_rate_raster_calculations.R](exceedance_rate_raster_calculations.R)
    * The calculations follow Davies et al. (2022) "From offshore to onshore PTHA via efficient Monte-Carlo sampling", for the case of stratified/importance sampling (Equations 17 and 20 in that paper).
    * An older code with similar functionality (but no variance calculations) is in [compute_exceedance_rates_for_threshold_depth_logic_tree_mean.R](compute_exceedance_rates_for_threshold_depth_logic_tree_mean.R). This has a finer-grained parallel approach which is much less efficient for the current problem, but in general it could be useful.
    * [compute_mean_exrate_upper_CI.R](compute_mean_exrate_upper_CI.R) Sums the logic-tree-mean exceedance-rate rasters on each source-zone, and also computes the upper limit of an approximate 95% confidence-interval for the all scenarios logic-tree-mean-exceedance-rate (NB: this does not account for epistemic uncertainties, just use of limited Monte-Carlo sampling). This makes various assumptions about the folder structure (I moved files to make things clean, see `make_directory_structure.sh`).

* [compute_exceedance_rates_for_threshold_stage_logic_tree_mean_newParallelPartition.R](compute_exceedance_rates_for_threshold_stage_logic_tree_mean_newParallelPartition.R) computes max-stage exceedance-rates, rather than depth exceedance-rates. 
    * It has an associated run script [run_compute_exceedance_rates_for_threshold_stage_logic_tree_mean_newParallelPartition.sh](run_compute_exceedance_rates_for_threshold_stage_logic_tree_mean_newParallelPartition.sh). 
    * After this script is run, we move the outputs to another directory with the script `move_max_stage_exceedance_rate_rasters_into_subfolder.sh`. Further post-processing scripts with documentation are found in the [directory where we move the outputs](ptha18-BunburyBusseltonRevised-sealevel60cm/max_stage_exceedance_rates/).

* [compute_exceedance_rates_at_epistemic_uncertainty_percentile.R](compute_exceedance_rates_at_epistemic_uncertainty_percentile.R) This computes the depth exceedance-rates at a given epistemic uncertainty percentile (e.g. 84%, or 16%). It only operates on one source-zone at a time, so results of source-zones must be later summed.
    * These percentile-uncertainty calculations are substantially more computationally demanding than the logic-tree-mean calculations. To spread the work over multiple PBS jobs, the script [make_exceedance_rate_jobs.R](make_exceedance_rate_jobs.R) makes a set of PBS scripts that collectively do all the work. It uses the template [run_compute_exceedance_rates_at_epistemic_uncertainty_RUNDIR_PERCENTILE_LOWER_UPPER.sh](run_compute_exceedance_rates_at_epistemic_uncertainty_RUNDIR_PERCENTILE_LOWER_UPPER.sh).
    * To sum the results over source-zones use [compute_sum_of_percentiles.R](compute_sum_of_percentiles.R). This makes various assumptions about the folder structure (in practice I moved tifs computed above to keep the structure clean, see `make_directory_structure.sh`).

* Scripts that begin with `make_vrt_` are used to create vrt files (which allow all domain tifs to be viewed in QGIS at once, as though they are a single file).

* [compute_exceedance_rates_for_multiple_threshold_depths_logic_tree_mean_newParallelPartition.R](compute_exceedance_rates_for_multiple_threshold_depths_logic_tree_mean_newParallelPartition.R) computes exceedance-rate rasters for a significant number of threshold depths.
    * It is run with [run_compute_exceedance_rates_for_multiple_threshold_depths_logic_tree_mean_newParallelPartition.sh](run_compute_exceedance_rates_for_multiple_threshold_depths_logic_tree_mean_newParallelPartition.sh)

## Calculation workflow

The workflow is not yet fully automated. The scripts include hard-coded assumptions, for example:
* The sources are sunda-arc thrust + outer-rise only. 
* The number of cores to use in parallel execution
* The name of the set of model runs (`BunburyBusseltonRevised`)

So the scripts will need modification for other applications and machines.

The order of running is basically like this...
```
#
# Do the logic-tree-mean and variance calculations for the thrust and outer-rise SundaArc sources on NCI
#
qsub run_compute_exceedance_rates_for_threshold_depth_logic_tree_mean_newParallelPartition.sh

#
# Do the epistemic uncertainty calculations for thrust and outer-rise sundaArc sources on NCI
#

# First make the job scripts (so we can spread calculations over nodes)
Rscript make_exceedance_rate_jobs.R

# Then submit all the resulting scripts
for i in run_compute_exceedance_rates_at_epistemic_uncertainty_ptha18-BunburyBusseltonRevised-sealevel60cm_*; do echo $i; qsub $i; done


#
# ... wait until the above runs have finished ...
#

# Once everything is complete, move the results to a local folder
source make_directory_structure.sh

#
# Combine thrust/outer-rise results for the logic-tree mean and Monte-Carlo variance.
#
cd ptha18-BunburyBusseltonRevised-sealevel60cm/highres_with_variance
Rscript ../../compute_mean_exrate_upper_CI.R
cd ptha18-BunburyBusseltonRevised-sealevel60cm-sum_of_source_zones/
source ../../../make_vrt_mean_variance_CI_summed_sources.sh
#
# We have now computed the logic-tree-mean exceedance-rate for 1mm inundation depth (sum of thrust and outer-rise) 
# and the upper limit of a 95% confidence interval for the latter (accounts for monte carlo errors due to limited sampling).
# They are in the folder where the last script was executed, see:
# - [Best estimate] raster summed_HS_exrate_all.vrt. 
# - [Upper limit of 95% Monte Carlo confidence interval] summed_HS_exrate_upperMonteCarloCI_all.vrt
#
# Go back to the main directory
cd ../../../

#
# Compute the sum of the epistemic uncertainty percentiles.
#
# This corresponds to the "co-monotonic" uncertainty treatment used in PTHA18.
# We deliberately skip domain 1 (global domain) because it requires too much
# memory -- to avoid this next-time I should not merge the global domain, but
# rather process it as separate tiles.
#
Rscript compute_sum_of_percentiles.R 84
Rscript compute_sum_of_percentiles.R 16

# ... wait until the above runs have finished ...
# Make vrts with the results, 16th percentile
cd ptha18-BunburyBusseltonRevised-sealevel60cm/highres_epistemic_uncertainty/16pc/ptha18-BunburyBusseltonRevised-sealevel60cm-depth_exrate_0.001_0.16_sum_of_source_zones/
source ../../../../make_vrt_percentiles.sh
# The 16th percentile result is in the folder where we executed the above script. 
# Its name is exceedance_rate_percentile_16pc.vrt
# It gives the 16th percentile of the epistemic uncertainty distribution for
# exceedance-rates of 1mm depth from the sunda arc thrust + outer-rise
# source-zones.
cd ../../../../

# Make vrts with the results, 84th percentile
cd ptha18-BunburyBusseltonRevised-sealevel60cm/highres_epistemic_uncertainty/84pc/ptha18-BunburyBusseltonRevised-sealevel60cm-depth_exrate_0.001_0.84_sum_of_source_zones/
source ../../../../make_vrt_percentiles.sh
# The 84th percentile result is in the folder where we executed the above script. 
# Its name is exceedance_rate_percentile_84pc.vrt
# It gives the 84th percentile of the epistemic uncertainty distribution for
# exceedance-rates of 1mm depth from the sunda arc thrust + outer-rise
# source-zones.
cd ../../../../


#
# Logic-tree mean calculations in terms of stage (rather than depth)
#
qsub run_compute_exceedance_rates_for_threshold_stage_logic_tree_mean_newParallelPartition.sh
# ... wait until the above run has finished, then do ...
source move_max_stage_exceedance_rate_rasters_into_subfolder.sh
# Now move to the folder ptha18-BunburyBusseltonRevised-sealevel60cm/max_stage_exceedance_rates/
# and post-process the results as described in the README therein.

```

## Tests

The script [exceedance_rate_raster_calculations.R](exceedance_rate_raster_calculations.R) has a function `.test_exceedance_rate_raster_calculations` which compares the Monte Carlo exceedance-rate and error variance  values computed with the main workhorse function herein, VS values computed independently using the model gauge outputs. This is for max-stage values of 1m above the background sea-level, at a particular site. I ran it and it works. Note these tests include hard-coded links to our specific simulations (so need to be updated if applying to another case).

The script [test_compute_exceedance_rates_at_epistemic_uncertainty.R](test_compute_exceedance_rates_at_epistemic_uncertainty.R) implements a test similar to the one above, for the case of exceedance-rates at an epistemic uncertainty percentile. I ran it and it works.

Another test that is useful is to compute a rasterised version of the max-stage exceedance-rates, and compare with PTHA18. These won't be the same, but should be similar in deep water far from the coast [and indeed they are].

Separately we can compare the Monte Carlo results at individual gauges with PTHA18; see calculations in [../max_stage_at_point](../max_stage_at_point). On theoretical grounds we expect them to agree well at sites far offshore where the linearity assumption works well. They are expected to disagree at sites in shallower water or closer to the coast, where the PTHA18 linear solver does not represent the hydrodynamics so well. In practice both of these theoretical expectations are met, with the former providing a useful check that the caculations have been correctly implemented.


