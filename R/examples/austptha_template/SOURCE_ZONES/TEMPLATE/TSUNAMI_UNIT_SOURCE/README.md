Run tsunami models for all unit-sources
---------------------------------------

The codes in this folder are made to run tsunami propagation models for every
unit-source on the source-zone. The unit-sources themselves must already have
been defined using code in [../EQ_SOURCE/](../EQ_SOURCE/).

**Currently the folder is set-up to run the SWALS model on raijin.nci.org.au. The
codes would require significant modification to be run on another machine, or using
another model.**

[create_files_for_all_simulations.R](create_files_for_all_simulations.R) is
used to make directory structures for all model runs (using some definitions
in [config.R](config.R)]. We make directories for all model runs by modifying 
template files in the [template](template) folder. This results in a separate
SWALS model that runs for each unit-source. 

[run_16.sh](run_16.sh) is used to submit 16 of the SWALS models (which are not
already slated for running) to the job queue, to be run on a single node. The
number 16 was chosen because that is efficient for the hardware on
raijin.nci.org.au.  Note that significant performance benefits were gained by
using 'numactl' when running multiple jobs on a single node. 

[run_permute_netcdf_gauges.sh](run_permute_netcdf_gauges.sh). This script
permutes the dimension of the tide-gauge netcdf output files, to ensure fast
access to single-station data [which is extremely convenient for later
analysis]. 


[check_runs_complete.R](check_runs_complete.R) and [check_logfile.R](check_logfile.R)
are used to check that the SWALS model runs have finished, and to do some basic
graphical QC checks. 
