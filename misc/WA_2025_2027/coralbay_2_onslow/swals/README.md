# Coral Bay to Onslow tsunami model

The steps to run are outlined below.

## 0. Make the model 
* Use the project modules from `source ../modules_ifx_2025_llvm.sh` to load modules needed to compile SWALS
* There are a few variants of the model to be compiled
  * Regular version: `make build`
  * Debug version - if models go unstable then this one will provide more information on the problematic location: `make debug`
  * Debug version with alternative nesting: `make both_debug_and_old_nesting_target`

## 1. Initialise and inspect the model
* Run a short test case (e.g. `qsub test/small_TEST_ONLY_initialise.pbs`) then once done:
  * Use `post-process/load_balance_script.R` to make a load balance file and inspect the elevation0 and max_stage rasters (or use ./post-process/process_initial_run.pbs to do both at once).
* Use [./post-process/report_domain_runtimes.R](./post-process/report_domain_runtimes.R) to inspect the domains with the worst run times. If one or a few domains take significantly longer than the rest consider coarsening them or splitting them. To split a domain, append 3 more or 8 more integers to it's line in the load balance file (to make 4 or 9 splits in the domain's rectangle). Then re-run SWALs using that load balance file. Then run the `load_balance_script.R` again on that output to load balance the partitions of this domain.
* Update the namelist files in [multidomain design control](./multidomain_design_control/) to point to this load balance file (find and replace does them all at once).
* Use [./post-process/make_rasters.R](./post-process/make_rasters) to make the elevation rasters. Inspect them in GIS, the bilinear interpolation may look different from what you expected. Adjust if required and re make the model.
* You can confirm the load balance wouldn't be improved by running `qsub test/small_TEST_ONLY_check_load_balance.pbs` and computing the load balance for this run. It should not have be able to improve significantly.
* It's also worth inspecting the tidally adjust elevation that the model will actually use. Pay attention to low lying areas that are disconnected from the ocean. You can patch these dry by setting the initial stage low.

## 2. Run the extreme and small test cases
* Qsub scripts in [./test](./test) from within the current (swals) directory
  * `qsub test/extreme.pbs`
  * `qsub test/extreme2.pbs`
  * `qsub test/small.pbs`
  * `qsub test/instability_carnarvon.pbs`

* After they have completed
  * Check mass conservation/energy decay with [./post-process/check_log_files.R](./post-process/check_log_files.R)
  * If sites exhibit velocity growth over time, investigate using the rasters (next).
  * Use [./post-process/make_rasters.R](./post-process/make_rasters.R) to make rasters with max-quantities of interest (as well as the last timestep UH or VH), then visualise in GIS and check for anything unexpected.
  * Use make rasters to inspect the uh and vh at the last timestep for the small case. It can reveal spurious flows from domain nesting nuances at nesting interfaces.

* Do convergence testing of the spatial resolution by submitting:
  * `qsub test/convergence_coarser_sumatra2004_fujisatake07.pbs` and
  * `qsub test/convergence_base_sumatra2004_fujisatake07.pbs`
  Inspect that they give similar solutions and timeseries.
* For a harder convergence test, look at the onshore results of a scenario with moderate inundation:
  * `qsub test/convergence_coarse_large_scenario.pbs`
  * `qsub test/convergence_base_large_scenario.pbs`
  * `qsub test/convergence_finer_large_scenario.pbs`

Proceed if you're satisfied with finalising the model's elevation and load balance.

## 3. Run the historic test cases
* Qsub scripts in [./historic](./historic)
  * `qsub historic/sumatra2004_fujisatake07.pbs`
  * `qsub historic/sumatra2004_fujii21.pbs` (optional)
  * `qsub historic/sumatra2005_Fujii20.pbs`
  * `qsub historic/sumatra2007_FujiiSatake08.pbs`
  * `qsub historic/java2006_fujiisatake2006_1kms_time_varying_with_horiz.pbs`

* After they have completed
  * Check mass conservation/energy decay with [./post-process/check_log_files.R](./post-process/check_log_files.R)
  * Check tide gauges with the `process_gauges_XXX.R` and `plot_gauges_XXX.R` scripts in the [../analysis/gauge_plots](../analysis/gauge_plots) directory. From inside that directory, create PNG files with model-vs-observed at tide gauges as follows: (here assuming the Java 2006 scenario)
    * `Rscript process_gauges_java2006.R ../swals/OUTPUTS/path-to-model-multidomain-dir-for-java2006` which makes an RDS file (`gauges_plot_XXXXXX.RDS`)
    * `Rscript plot_gauges_java2006.R ../swals/OUTPUTS/path-to-model-multidomain-dir-for-java2006/gauges_plot_XXXXX.RDS 0.5 1 6` which plots the model and observations, here with y-axes covering +- 0.5m and x-axes (time) covering 1 to 6 hours post-earthquake.
  * Use [./post-process/make_rasters.R](./post-process/make_rasters.R) to make rasters with max-quantities of interest, then visualise in GIS and check for anything unexpected

## 4. Create the PBS scripts to run relevant scenarios in [../sources/](../sources/), and then submit them
Use `Rscript create_random_ptha_qsub_scripts_hazard.R` to make scripts to run the hazard scenarios.

## 5. After the submitted scripts have run
Check they finished and that the mass/energy/etc are OK
  * This can partly use [./post-process/check_log_files.R](./post-process/check_log_files.R), but also make a note of how many runs should have completed (i.e. how many sources were to be modelled), since the latter script won't show if a run is completely missing.

## 6. Debug any failing runs.
* If a run failed, then its log files may be useful for understanding the cause of the error
  * Sometimes the issues are hardware failures, so the model can be re-run as-is without problems
  * A few scenarios failed due to intense inundation around NW Cape, but could be made to run stably by reducing the timestep. 
    * This can be achieved using [multidomain_design_control/hazard_lowts.nml](./multidomain_design_control/low_ts.nml) (instead of [hazard.nml](./multidomain_design_control/hazard.nml)) as the configuration file. 
    * Although a single timestep reduction is hard-coded in the latter file, in practice we tried a sequence of timestep reductions (3, then 5 for those that still failed, then 7). A number of initially failing scenarios worked using a 3x timestep reduction. For those that still failed, some worked with 5x. One needed 7x. 
* Tips for rerunning models
  * Before re-running the model, you must delete its multidomain directory
  * After the new runs complete, repeat step 4 to check all the runs are OK, and if needed continue debugging.
* Scripts to make qsub scripts for the particular models that failed for me are in [post-process/make_jobs_failed_runs_using_lowts.R](post-process/make_jobs_failed_runs_using_lowts.R) and [post-process/make_jobs_failed_runs_using_alternative_nesting](post-process/make_jobs_failed_runs_using_alternative_nesting). 

## 7. Create rasters
From inside the `post-process` directory and submit `qsub create_tarred_rasters_from_tarred_multidomains.pbs`

## 8. Compress the model output tar files (basename like `RUN_(DATETIME_STAMP).tar`).
I did this by editing `post-process/bzip2_some_files.R` and running it with the PBS run script in the same directory.

## Done
Now go to the [../analysis](../analysis/) folder.
