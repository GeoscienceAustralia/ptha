# Kalbarri to Coral Bay tsunami model

The steps to run are outlined below.

0. Make the model 
* Use `source R_431_NCI_modules.sh` to load modules needed to compile SWALS
* There are a few variants of the model to be compiled
  * Regular version: `make build`
  * Debug version - if models go unstable then this one will provide more information on the problematic location: `make debug`
  * Debug version with alternative nesting: `make both_debug_and_old_nesting_target`
* Run a short test case (e.g. `qsub run_tests/run_smallTestCase_TESTONLY.sh`) and use the results to make a load balance file, using `post_process/load_balance_script.R`.

1. Run the realistic test cases
* Qsub scripts in [./run_tests](./run_tests)
  * `qsub run_tests/run_sumatra2004_kalbarri2coralbay_fujisatake07.sh`
  * `qsub run_tests/run_sumatra2004_kalbarri2coralbay_fujii21.sh`
  * `qsub run_tests/run_sumatra2005_kalbarri2coralbay_Fujii20.sh`
  * `qsub run_tests/run_sumatra2007_kalbarri2coralbay_FujiiSatake08.sh`
  * `qsub run_tests/run_java2006_kalbarri2coralbay_fujiisatake2006_1kms_time_varying_with_horiz.sh`
* After they have completed
  * Check mass conservation/energy decay with [./post_process/check_log_files.R](./post_process/check_log_files.R)
  * Check tide gauges with the `process_gauges_XXX.R` and `plot_gauges_XXX.R` scripts in the [./plots](./plots) directory. From inside that directory, use comments like below to create PNG files with model-vs-observed at tide gauges (here assuming the Java 2006 scenario)
    * `Rscript process_gauges_java2006.R ../OUTPUTS/path-to-model-multidomain-dir-for-java2006` which makes an RDS file (`gauges_plot_XXXXXX.RDS`)
    * `Rscript plot_gauges_java2006.R ../OUTPUTS/path-to-model-multidomain-dir-for-java2006/gauges_plot_XXXXX.RDS 0.5 1 6` which plots the model and observations, here with y-axes covering +- 0.5m and x-axes (time) covering 1 to 6 hours post-earthquake.
  * Use [./post_process/make_rasters.R](./post_process/make_rasters) to make rasters with max-quantities of interest, then visualise in GIS and check for anything unexpected

2. Run the extreme and small test cases
* Qsub scripts in [./run_scripts](./run_scripts) from within the current directory
  * `qsub run_scripts/run_extremeTestCase.sh`
  * `qsub run_scripts/run_extremeTestCase2.sh`
  * `qsub run_scripts/run_smallTestCase.sh`
* After they have completed
  * Check mass conservation/energy decay with [./post_process/check_log_files.R](./post_process/check_log_files.R)  
    * These show a few sites with velocity growth over time that were investigated using the rasters (below).
  * Use [./post_process/make_rasters.R](./post_process/make_rasters) to make rasters with max-quantities of interest (as well as the last timestep UH or VH), then visualise in GIS and check for anything unexpected.
    * These were used to check a few sites with velocity growth over time. They occur in very localised spots at a couple of nesting boundaries, but did not appear to be having an important effect. 

3. Create the PBS scripts to run relevant scenarios in [../sources/](../sources/), and then submit them
  * Use `Rscript create_random_ptha_qsub_scripts_hazard.R` to make scripts to run the hazard scenarios. 

4. After the submitted scripts have run
  * Check they finished and that the mass/energy/etc are OK
    * This can partly use [./post_process/check_log_files.R](./post_process/check_log_files.R), but also make a note of how many runs should have completed (i.e. how many sources were to be modelled), since the latter script won't show if a run is completely missing.

5. Debug any failing runs.
  * If a run failed, then its log files may be useful for understanding the cause of the error
    * Sometimes the issues are hardware failures, so the model can be re-run as-is without problems
    * A few scenarios failed due to intense inundation around NW Cape, but could be made to run stably by reducing the timestep. 
      * This can be achieved using multidomain_kalbarri2onslow_20251218_hazard_lowts.nml (instead of multidomain_kalbarri2onslow_20251218_hazard.nml) as the configuration file. 
      * Although a single timestep reduction is hard-coded in the latter file, in practice we tried a sequence of timestep reductions (3, then 5 for those that still failed, then 7). A number of initially failing scenarios worked using a 3x timestep reduction. For those that still failed, some worked with 5x. One needed 7x. 
  * Tips for rerunning models
    * Before re-running the model, you must delete its multidomain directory
    * After the new runs complete, repeat step 4 to check all the runs are OK, and if needed continue debugging.
  * Scripts to make qsub scripts for the particular models that failed for me are in [post_process/make_jobs_failed_runs_using_lowts.R](post_process/make_jobs_failed_runs_using_lowts.R) and [post_process/make_jobs_failed_runs_using_alternative_nesting](post_process/make_jobs_failed_runs_using_alternative_nesting). 

6. Create rasters by going inside the `post_process` directory and doing `qsub create_tarred_rasters_from_tarred_multidomains.pbs`

7. Compress the model output tar files (basename like `RUN_(DATETIME_STAMP).tar`).
  * I did this by editing `post_process/bzip2_some_files.R` and running it with the PBS run script in the same directory.
