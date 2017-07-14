Code for making earthquake and tsunami events
---------------------------------------------

Codes in this folder should be run after all unit-source tsunami runs are
complete (see [../TSUNAMI_UNIT_SOURCES](tsunami_unit_sources) ).

The code is designed to work on raijin.nci.org.au, using the SWALS model, 
and will need to be ported to work on any other platform / or with any other 
solver. 


# Step 1

Define key parameters by modifying [config.R](config.R). These parameters
control e.g. the range of earthquake magnitudes that our tsunami events cover, 
the number of stochastic slip events created, computational options, among
other things. 

# Step 2

To make earthquake and tsunami events (for both stochastic and uniform slip),
see [run_make_all_tsunami_events.PBS](run_make_all_tsunami_events.PBS). 

    qsub run_make_all_tsunami_events.PBS

# Step 3

On large source-zones, the final stochastic_slip generation step (from Step 2)
might not complete [say if the above script is killed due to exceeding time
limits]. If that happens, you need to follow Step 3 (otherwise skip to Step 4). 

The way to tell if this has happened is to examine raijin's job files (.o and
.e files generated when the job ends). If the job was killed, they will
contain statements to that effect. 

If the job was killed, then not all of your stochastic slip tsunami have been created. 
You can create the unfinished stochastic slip with
[run_unfinished_stochastic_earthquake_tsunami.R](run_unfinished_stochastic_earthquake_tsunami.R).

    Rscript run_unfinished_stochastic_earthquake_tsunami.R

This will split the remaining work over multiple jobs (so it is reasonably
fast), and submit all jobs to the PBS queue. Each of those jobs will make an
RDS file with its portion of the unfinished results. Those are not initially
written to the main netCDF output file, to avoid the chance of multiple files
trying to write at once (which is unsupported and can lead to data loss).
However, once all the 'unfinished' runs are finished, we put the data into the
main netcdf file using the script
[run_merge_additional_stochastic_tsunami_into_netcdf.PBS](run_merge_additional_stochastic_tsunami_into_netcdf.PBS)

    qsub run_merge_additional_stochastic_tsunami_into_netcdf.PBS

Once that has finished and you have confirmed the netCDF files are ok, you
might want to delete the folder R_images_tmp that was used to store temporary
RDS files. 

