R scripts for processing the inputs and outputs for SWALS.

* [create_plots_from_tarred_multidomain_dirs.R](create_plots_from_tarred_multidomain_dirs.R) can make basic png images of various flow maxima (stage, speed, flux) from the tarred multidomain directories.
* [create_tarred_rasters_from_tarred_multidomains.R](create_tarred_rasters_from_tarred_multidomains.R) makes rasters for an existing tarred multidomain folder, saving them to another tar archive. This approach of using tarfiles is important to avoid making too many files (violating our iinode quota on NCI).
* [make_rasters.R](make_rasters.R) and [make_all_rasters.sh](make_all_rasters.sh) are useful for making rasters from a single model run.
* [load_balance_script.R](load_balance_script.R) can be run from inside a multidomain directory, and will produce a file `load_balance_file.txt` which tries to distribute domains among MPI images to equally distribute the work from that run.
    * The text files in [../multidomain_design_control](../multidomain_design_control) were produced this way (but apply to different model setups and core-counts).
* [make_domains_shapefile.R](make_domains_shapefile.R) can make a shapefile depicting the multidomain layout.
* [report_domain_runtimes.R](report_domain_runtimes.R) can summarise the time required by each domain, which is useful for optimization.
* The script [check_log_files.R](check_log_files.R) runs a number of basic QC checks on the tsunami inundation simulations by scanning the log file. See notes in the script header for details of the theory.
