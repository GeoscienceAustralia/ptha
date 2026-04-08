# Exploratory model covering the full Kalbarri to Onslow region

This is not as up-to-date as some subsequent models in the region, but was used to de-risk the analysis and to test the tidal adjustment technique.


1. Make the model with `make build`

2. Run the realistic test cases
* Qsub scripts in [./run_scripts](./run_scripts) 
   * **(FIXME clean this up, there are many experiments here too. Also, as of Jan 2026 the validation tests were not made for the most up-to-date models, I ran them with earlier versions (will be fine, but update))**
* After they have completed
  * Check mass conservation/energy decay with [./post_process/check_log_files.R](./post_process/check_log_files.R)
  * Check tide gauges with the `process_gauges_XXX.R` and `plot_gauges_XXX.R` scripts in the [./plots](./plots) directory. From inside that directory, use comments like below to create PNG files with model-vs-observed at tide gauges (here assuming the Java 2006 scenario)
    * `Rscript process_gauges_java2006.R ../OUTPUTS/path-to-model-multidomain-dir-for-java2006` which makes an RDS file (`gauges_plot_XXXXXX.RDS`)
    * `Rscript plot_gauges_java2006.R ../OUTPUTS/path-to-model-multidomain-dir-for-java2006/gauges_plot_XXXXX.RDS 0.5 1 6` which plots the model and observations, here with y-axes covering +- 0.5m and x-axes (time) covering 1 to 6 hours post-earthquake.
  * Use [./post_process/make_rasters.R](./post_process/make_rasters) to make rasters with max-quantities of interest, then visualise in GIS and check for anything unexpected

3. Run the extreme test cases
* Qsub scripts in [./run_scripts](./run_scripts) from within the current directory **(FIXME clean this up, be specific)**
  * `qsub run_scripts/run_extremeTestCase_kalbarri2onslow_20251218.sh`
  * `qsub run_scripts/run_smallTestCase_kalbarri2onslow_20251218.sh`
* After they have completed
  * Check mass conservation/energy decay with [./post_process/check_log_files.R](./post_process/check_log_files.R)  
    * **(FIXME: The `dry-runs` had a later time growth in max-speed that warrents investigation before doing the final runs)**
  * Use [./post_process/make_rasters.R](./post_process/make_rasters) to make rasters with max-quantities of interest (as well as the last timestep UH or VH), then visualise in GIS and check for anything unexpected.

3. Create the PBS scripts to run relevant scenarios in [../sources/](../sources/), and then submit them
  * e.g. `Rscript create_random_ptha_qsub_scripts_hazard.R` to make scripts to run the hazard scenarios. 
  * An alternative set of scenarios which include the effect of sloping bathymetry (horizontal components) in the PTHA18 source models are made with `Rscript create_random_ptha_qsub_scripts_hazard_bathyslope.R`. 

4. After the submitted scripts have run
  * Check they finished and that the mass/energy/etc are OK
    * This can partly use [./post_process/check_log_files.R](./post_process/check_log_files.R), but also make a note of how many runs should have completed (i.e. how many sources were to be modelled), since the latter script won't show if a run is completely missing.

5. Debug any failing runs.
  * If a run failed, then its log files may be useful for understanding the cause of the error
    * Sometimes the issues are hardware failures, so the model can be re-run as-is without problems
    * A few scenarios failed due to intense inundation around NW Cape, but could be made to run stably by reducing the timestep. 
      * This can be achieved using multidomain_kalbarri2onslow_20251218_hazard_lowts.nml (instead of multidomain_kalbarri2onslow_20251218_hazard.nml) as the configuration file. 
      * Although a single timestep reduction is hard-coded in the latter file, in practice we tried a sequence of timestep reductions (3, then 5 for those that still failed, then 7). A number of initially failing scenarios worked using a 3x timestep reduction. For those that still failed, some worked with 5x. One needed 7x. 
  * Before re-running the model, you must delete its multidomain directory
  * After the new runs complete, repeat step 4 to check all the runs are OK, and if needed continue debugging.

6. Create rasters by going inside the `post_process` directory and doing `qsub create_tarred_rasters_from_tarred_multidomains.pbs`

7. Compress the model output tar files (basename like `RUN_(DATETIME_STAMP).tar`) 
