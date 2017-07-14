Run tsunami models for all unit-sources
---------------------------------------

The codes in this folder are made to run tsunami propagation models for every
unit-source on the source-zone. The unit-sources themselves must already have
been defined using code in [../EQ_SOURCE/](../EQ_SOURCE/).

**Currently the folder is set-up to run the SWALS model on raijin.nci.org.au. The
codes would require significant modification to be run on another machine, or using
another model.**

[create_files_for_all_simulations.R](create_files_for_all_simulations.R) is
used to make directory structures for all model runs. This modifies files in
the [template](template) folder, to make a SWALS model that runs for each
unit-source. Because it is modifying SWALS input files, it is obviously strongly
dependent on the input-file design for the SWALS model. Furthermore, to run
the code on raijin.nci.org.au we need to load particular modules and set
particular paths. This code takes care of that, but it would require modification
to work on systems other than raijin.nci.org.au.

[run_16.sh](run_16.sh) is used to submit 16 of the SWALS models (which are not
already slated for running) to the job queue. This is set-up to work on
raijin.nci.org.au. Note that significant performance benefits were gained by using
'numactl' when running multiple jobs on a single node. 

[run_permute_netcdf_gauges.sh](run_permute_netcdf_gauges.sh). This script
permutes the dimension of the tide-gauge netcdf output files, to ensure fast
access to single-station data [which is extremely convenient for later
analysis]. It is set-up for raijin.nci.org.au


[check_runs_complete.R](check_runs_complete.R) and [check_logfile.R](check_logfile.R)
are used to check that the SWALS model runs have finished, and to do some basic
graphical QC checks. 
