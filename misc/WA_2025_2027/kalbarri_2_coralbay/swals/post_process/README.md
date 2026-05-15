R scripts for processing the inputs and outputs for SWALS.

They contain a hard-coded reference to the `ptha` repository. To change it quickly, see [HOW_TO_QUICKLY_CHANGE_SWALS_PLOT_CODE_PATH.txt](HOW_TO_QUICKLY_CHANGE_SWALS_PLOT_CODE_PATH.txt).

* [check_log_files.R](check_log_files.R) is used to check that runs completed and have the expected mass and energy conservation.
* [create_plots_from_tarred_multidomain_dirs.R](create_plots_from_tarred_multidomain_dirs.R) can make basic png images of various flow maxima (stage, speed, flux) from the tarred multidomain directories.
* [create_tarred_rasters_from_tarred_multidomains.R](create_tarred_rasters_from_tarred_multidomains.R) makes rasters for an existing tarred multidomain folder, saving them to another tar archive. This approach of using tarfiles is important to avoid making too many files (violating our iinode quota on NCI).
* [make_rasters.R](make_rasters.R) and [make_all_rasters.sh](make_all_rasters.sh) are useful for making rasters from a single model run.
* Scripts to run jobs that failed the first time. After running all the tsunami scenarios, I use [check_log_files.R](check_log_files.R) to help find runs that didn't finish or might otherwise be problematic. To fix them, they were run with one or more alternative setups until they completed successfully.
  * [make_jobs_failed_runs_using_alternate_nesting.R](make_jobs_failed_runs_using_alternate_nesting.R) is the first alternative approach that is used. This runs the models using a slightly different nesting scheme, in `debug` mode which means that if they fail then I will get better information on where the failures occurred.
  * [make_jobs_failed_runs_using_lowts.R](make_jobs_failed_runs_using_lowts.R) was used to run models that didn't work under the previous approach, using a reduced (1/3) timestep.  
  * One job still failed, but that one worked with an even more reduced timestep (1/7) -- which I setup manually. 
* [load_balance_script.R](load_balance_script.R) can be run from inside a multidomain directory, and will produce a file `load_balance_file.txt` which tries to distribute domains among MPI images to equally distribute the work from that run.
    * The text file in [../multidomain_design_control](../multidomain_design_control) was produced this way.
* [make_domains_shapefile.R](make_domains_shapefile.R) can make a shapefile depicting the multidomain layout.
* [report_domain_runtimes.R](report_domain_runtimes.R) can summarise the time required by each domain. This is useful to identify slow domains (for optimization).
* The script [check_log_files.R](check_log_files.R) runs a number of basic QC checks on the tsunami inundation simulations by scanning the log file. See notes in the script header for details of the theory.
* [run_bzip2_some_files.sh](run_bzip2_some_files.sh) will run the script [bzip2_some_files.R](bzip2_some_files.R) -- make sure you edit the latter to point to the correct files.
