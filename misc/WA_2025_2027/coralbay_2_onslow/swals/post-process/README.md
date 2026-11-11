R scripts for processing the inputs and outputs for SWALS.

They contain a hard-coded reference to the `ptha` repository.

* [check_log_files.R](check_log_files.R) is used to check that runs completed and have the expected mass and energy conservation.
* [create_tarred_rasters_from_tarred_multidomains.R](create_tarred_rasters_from_tarred_multidomains.R) makes rasters for an existing tarred multidomain folder, saving them to another tar archive. This approach of using tarfiles is important to avoid making too many files (violating our iinode quota on NCI).
* [make_rasters.R](make_rasters.R) is useful for making rasters from a single model run.
* Scripts to run jobs that failed the first time. After running all the tsunami scenarios, I use [check_log_files.R](check_log_files.R) to help find runs that didn't finish or might otherwise be problematic. To fix them, they were run with one or more alternative setups until they completed successfully.
* [make_jobs_failed_runs_using_lowts.R](make_jobs_failed_runs_using_lowts.R) was used to run models that didn't work under the previous approach, using a reduced (1/3) timestep.
* Jobs still failing were rerun with an even more reduced timestep (1/7) using [make_jobs_failed_runs_using_lowerts.R](make_jobs_failed_runs_using_lowerts.R).
* [load_balance_script.R](load_balance_script.R) can be run from inside a multidomain directory, and will produce a file `load_balance_file.txt` which tries to distribute domains among MPI images to equally distribute the work from that run.
    * The text file in [../multidomain_design_control](../multidomain_design_control) was produced this way.
* [make_domains_shapefile.R](make_domains_shapefile.R) can make a shapefile depicting the multidomain layout.
* [report_domain_runtimes.R](report_domain_runtimes.R) can summarise the time required by each domain. This is useful to identify slow domains (for optimization).
* The script [check_log_files.R](check_log_files.R) runs a number of basic QC checks on the tsunami inundation simulations by scanning the log file. See notes in the script header for details of the theory.
* [bzip2_some_files.pbs](bzip2_some_files.pbs) will run the script [bzip2_some_files.R](bzip2_some_files.R) -- make sure you edit the latter to point to the correct files.
* [process_initial_run.pbs](./process_initial_run.pbs) automates making the domain shapefiles, rasters, virtual rasters and load balancing for any simulation. Hard code the simulation's path in the script. This can be set to run after the simulation job finishes using pbs flags: `qsub -W depend=afterok:JOBID process_initial_run.pbs` where `JOBID` is the job id of the simulation job.
