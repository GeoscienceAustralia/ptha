Code for making earthquake and tsunami events
---------------------------------------------

Codes in this folder should be run after all unit-source tsunami runs are
complete (see [../TSUNAMI_UNIT_SOURCES](../TSUNAMI_UNIT_SOURCES) ).

The code is designed to work on raijin.nci.org.au, using the SWALS model, 
and will need to be ported to work on any other platform / or with any other 
solver. 


# Step 1

Define key parameters by modifying [config.R](config.R). These parameters
control e.g. the range of earthquake magnitudes that our tsunami events cover, 
the number of stochastic slip events created, computational options, among
other things. 

Also, ensure that the parameters in the user input table
[../../../DATA/SOURCEZONE_PARAMETERS/sourcezone_parameters.csv](../../../DATA/SOURCEZONE_PARAMETERS/sourcezone_parameters.csv)
are up to date (see the [README][../../../DATA/SOURCEZONE_PARAMETERS/README.md]
in that folder for more information).

# Step 2

To make earthquake and tsunami events (for both stochastic, variable_uniform,
and uniform slip), see
[run_make_all_tsunami_events.PBS](run_make_all_tsunami_events.PBS). 

    qsub run_make_all_tsunami_events.PBS

# Step 3

On large source-zones, the final stochastic_slip generation step (from Step 2)
might not complete [say if the above script is killed due to exceeding time
limits]. If that happens, you need to follow Step 3 - otherwise, skip now to Step 4. 

The way to tell if this has happened is to examine raijin's job files (.o and
.e files generated when the job ends). If the job was killed, they will
contain statements to that effect. 

If the job was killed, then not all of your stochastic slip tsunami have been
created.  You can create the unfinished stochastic slip with 
[run_unfinished_stochastic_earthquake_tsunami.R](run_unfinished_stochastic_earthquake_tsunami.R).
(Note: Currently this is automatically executed by Step 2, IF there are more than some
threshold number of events in stochastic event set):

    Rscript run_unfinished_stochastic_earthquake_tsunami.R --stochastic_slip

and similarly for the variable_uniform slip events (also will be automatically
executed in Step 2 if there are many events):

    Rscript run_unfinished_stochastic_earthquake_tsunami.R --variable_uniform_slip

The above commands split the remaining work over multiple PBS jobs (so it is reasonably
fast), and submit all jobs to the PBS queue. Each of those jobs will make an
RDS file with its portion of the unfinished results. Those are not initially
written to the main netCDF output file, to avoid the chance of multiple processes
trying to write to the same file at once (which is unsupported and can lead to data loss).

Once all the 'unfinished' runs are finished, we put the data into the
main netcdf file using the script
[run_merge_additional_stochastic_tsunami_into_netcdf.PBS](run_merge_additional_stochastic_tsunami_into_netcdf.PBS)

    qsub run_merge_additional_stochastic_tsunami_into_netcdf.PBS

This has to be executed manually. Once that has finished and you have confirmed
the netCDF files are ok, you might want to delete the folder R_images_tmp that
was used to store temporary RDS files. 

On large source-zones, in some instances it is possible for the unit-source, gauge
and earthquake event properties to not be correctly written to the files of the form
all_SLIP_TYPE_slip_earthquake_tsunami_SOURCEZONE.nc files. This happens if R
runs out of memory just after creating those files (i.e. so the file exists, but the
variables are not written to it). You might not even notice this until running
test code below. If this happens, then a 'quick fix' is to re-write the data to
the file separately (without affecting any flow variables). That can be done
with a command like the following (here for variable_uniform slip):

    Rscript make_all_earthquake_tsunami.R --variable_uniform_slip --only_update_non_flow_variables


# Step 4

At this point all the stochastic and uniform slip events have been created, and it is 
time to check the results. 

The script [check_events.R](check_events.R) can be used to plot properties of the earthquake
events, and sanity check peak wave heights, earthquake magnitudes, earthquake dimensions, etc.
Failures might indicate problems in writing the output files -- see above for possible fixes.

The script [check_dart_example.R](check_dart_example.R) gives an example of comparing
model scenarios to tide gauge data -- but obviously the file paths etc will need to be
adapted on a case-by-case basis. There are further codes in [plots](plots) which can
be used to examine the events. 

