Run tsunami models for all unit-sources
---------------------------------------

The codes in this folder are made to run tsunami propagation models for every
unit-source on the source-zone. The unit-sources themselves must already have
been defined using code in [../EQ_SOURCE/](../EQ_SOURCE/).

**Currently the folder is set-up to run the SWALS model on raijin.nci.org.au. The
codes would require significant modification to be run on another machine, or using
another model.**

To run the models, the process is

# Step 1:

Check the variables defined in [config.R](config.R) and see that they point
to the correct unit sources, and that the model output files will go to a location
that you are happy with.

# Step 2

[create_files_for_all_simulations.R](create_files_for_all_simulations.R) is
used to make directory structures for all model runs (using some definitions
in [config.R](config.R)). We make directories for all model runs by modifying 
template files in the [template](template) folder. This results in a separate
SWALS model that runs for each unit-source. It is run with:

    Rscript create_files_for_all_simulations.R

and is quite light-weight (i.e. might not need to use the job queue).

# Step 3

[run_16.sh](run_16.sh) is used to submit 16 of the SWALS models (which are not
already slated for running) to the job queue, to be run on a single node. The
number 16 was chosen because that is efficient for the hardware on
raijin.nci.org.au.  Note that significant performance benefits were gained by
using 'numactl' when running multiple jobs on a single node. It is run with:
    
    qsub run_16.sh

The above needs to be done repeatedly, until all model runs are completed.

# Step 4

Perform some basic checks that the models have run correctly.

[check_runs_complete.R](check_runs_complete.R) is used to check that the SWALS
model runs have finished, and to do some basic checks. It is run *from within
R* using: 
    
    source(check_runs_complete.R)
    # This is used to check for models that seem to have not completed
    check_models_have_been_run()
    # If the above is ok, then try this to check that gauge outputs exist and
    # are correctly ordered
    check_models_gauge_integrity()

[check_logfile.R](check_logfile.R) is used to plot the evolution of peak-stage
and mass balance in every tsunami model, using the log-files. To run it, do
    
    Rscript check_logfile.R

and then investigate the resulting pdf file 'tsunami_log_check.pdf'

# Step 5

[run_permute_netcdf_gauges.sh](run_permute_netcdf_gauges.sh) is used 
permutes the dimension of the tide-gauge netcdf output files, to ensure fast
access to single-station data [which is extremely convenient for later
analysis]. It is run with:

    qsub run_permute_netcdf_gauges.sh


