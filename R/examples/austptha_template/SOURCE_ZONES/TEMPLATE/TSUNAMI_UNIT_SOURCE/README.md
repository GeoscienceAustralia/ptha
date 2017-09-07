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
in [config.R](config.R)). 

We make directories for all model runs by modifying template files in the
[template](template) folder. This allows us to define separate output
directories and initial conditions (i.e. tsunami unit sources) for every model
run. Note that files in the [template](template) folder define other model
configuration options (e.g. resolution), which may need to be modified
depending on your applications.

Once the [template](template) files are ok, to create a separate SWALS model
that runs for each unit-source, do: 

    Rscript create_files_for_all_simulations.R

This script is quite light-weight (i.e. you might not need to use the job queue).

# Step 3

[run_16.sh](run_16.sh) is used to submit 16 of the SWALS models (which are not
already slated for running) to the job queue, to be run on a single node. The
number 16 was chosen because that is efficient for the hardware on
raijin.nci.org.au.  Note that significant performance benefits were gained by
using 'numactl' when running multiple jobs on a single node. It is run with:
    
    qsub run_16.sh

The above needs to be done repeatedly, until all model runs are completed.

**NOTE**: If you are running multiple source-zones, it may be preferable to use the
alternative script [../../run_16.sh](../../run_16.sh). The latter is run from
it's own directory, and has the advantage of being able to run models from more
than one source-zone (if it doesn't find 16 models to run on a single
source-zone).

# Step 4

Perform some basic checks that the models have run correctly.

[check_runs_complete.R](check_runs_complete.R) is used to check that the SWALS
model runs have finished, and to do some basic checks. It is run *from within
R* using: 
   
```r 
    # Read the functions into R
    source('check_runs_complete.R')
    # This is used to check for models that seem to have not completed
    check_models_have_been_run()
    # If the above is ok, then try this to check that gauge outputs exist and
    # are correctly ordered
    check_model_gauge_integrity()
```

Note that this process can be run over multiple source-zones at once, by using
the alternative script [../../checkruns.R](../../checkruns.R).


[check_logfile.R](check_logfile.R) is used to plot the evolution of peak-stage
and mass balance in every tsunami model, using the log-files. To run it, do
    
    Rscript check_logfile.R

and then visually investigate the resulting pdf file 'tsunami_log_check.pdf',
which will show the time-evolution of peak-stage and mass-balance, for every
model.

# Step 5

[run_permute_netcdf_gauges.sh](run_permute_netcdf_gauges.sh) is used 
permutes the dimension of the tide-gauge netcdf output files, to ensure fast
access to single-station data [which is extremely convenient for later
analysis]. It is run with:

    qsub run_permute_netcdf_gauges.sh


