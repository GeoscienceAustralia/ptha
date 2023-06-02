# Probabilistic Inundation calculations

An updated version of these codes exists in the "BunburyBusselton" folder. That version is better structured and better documented, so preferable for further work.

This folder contains scripts for probabilistic inundation calculation, which can be run after all the models in [../../swals](../../swals) have been simulated. 

## Key files

* [compute_exceedance_rates_for_threshold_depth_logic_tree_mean_newParallelPartition.R](compute_exceedance_rates_for_threshold_depth_logic_tree_mean_newParallelPartition.R) - This computes the logic-tree mean exceedance-rate for a specified depth threshold, as well as the variance of the Monte-Carlo error. 
    * run with [run_compute_exceedance_rates_for_threshold_depth_logic_tree_mean_newParallelPartition.sh](run_compute_exceedance_rates_for_threshold_depth_logic_tree_mean_newParallelPartition.sh)
      * This treats each source-zone separately. 
      * Results of source-zones must be later summed with [compute_mean_exrate_upper_CI.R](compute_mean_exrate_upper_CI.R). This also computes the upper limit of an approximate 95% confidence-interval for the all scenarios logic-tree-mean-exceedance-rate. The confidence interval does not account for epistemic uncertainties, just use of limited Monte-Carlo sampling.
    * The computational core underlying these calculations is in [exceedance_rate_raster_calculations.R](exceedance_rate_raster_calculations.R).
    * The calculations follow Davies et al. (2022) "From offshore to onshore PTHA via efficient Monte-Carlo sampling", for the case of stratified/importance sampling (Equations 17 and 20 in that paper).
    * An older code with similar functionality (but no variance calculations) is in [compute_exceedance_rates_for_threshold_depth_logic_tree_mean.R](compute_exceedance_rates_for_threshold_depth_logic_tree_mean.R). This has a finer-grained parallel approach which is much less efficient for the current problem, but in general could be useful.

* [compute_exceedance_rates_for_threshold_stage_logic_tree_mean_newParallelPartition.R](compute_exceedance_rates_for_threshold_stage_logic_tree_mean_newParallelPartition.R) is similar to the above, but computes max-stage exceedance-rates (rather than depth exceedance-rates). 
    * After calculation the results were manually moved to [reviseddomains_080422/exrates_for_multiple_stages/](reviseddomains_080422/exrates_for_multiple_stages/) for further processing, as documented therein.

* [compute_exceedance_rates_at_epistemic_uncertainty_percentile.R](compute_exceedance_rates_at_epistemic_uncertainty_percentile.R) This computes the depth exceedance-rates at a given epistemic uncertainty percentile (e.g. 84%, or 16%). It only operates on one source-zone at a time, so results of source-zones must be later summed.
    * These percentile-uncertainty calculations are substantially more computationally demanding than the logic-tree-mean calculations. To spread the work over multiple PBS jobs, the script [make_exceedance_rate_jobs.R](make_exceedance_rate_jobs.R) makes a set of PBS scripts that collectively do all the work (BEWARE: Thus currently uses a hard-coded number of domains, it should be edited to determine this automatically). It uses the template [run_compute_exceedance_rates_at_epistemic_uncertainty_RUNDIR_PERCENTILE_LOWER_UPPER.sh](run_compute_exceedance_rates_at_epistemic_uncertainty_RUNDIR_PERCENTILE_LOWER_UPPER.sh).
    * To sum the results over source-zones use [compute_sum_of_percentiles.R](compute_sum_of_percentiles.R). 

* Scripts that begin with `make_vrt_` are used to create vrt files (which allow all domain tifs to be viewed in QGIS at once, as though they are a single file).

* [compute_exceedance_rates_for_multiple_threshold_depths_logic_tree_mean_newParallelPartition.R](compute_exceedance_rates_for_multiple_threshold_depths_logic_tree_mean_newParallelPartition.R) computes exceedance-rate rasters for a significant number of threshold depths.
    * It is run with [run_compute_exceedance_rates_for_multiple_threshold_depths_logic_tree_mean_newParallelPartition.sh](run_compute_exceedance_rates_for_multiple_threshold_depths_logic_tree_mean_newParallelPartition.sh)
    * After this is completed, the resulting folders were manually moved to [reviseddomains_080422/exrates_for_multiple_depths/](reviseddomains_080422/exrates_for_multiple_depths/) for further processing. See the README therein for details.
