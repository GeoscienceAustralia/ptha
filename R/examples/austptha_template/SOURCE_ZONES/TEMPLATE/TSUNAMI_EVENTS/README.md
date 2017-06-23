Code for making earthquake and tsunami events
---------------------------------------------

Codes in this folder should be run after all unit-source tsunami runs are
complete (see [../TSUNAMI_UNIT_SOURCES](tsunami_unit_sources) ).

To make earthquake and tsunami events (for both stochastic and uniform slip),
see [run_make_all_tsunami_events.PBS](run_make_all_tsunami_events.PBS). 

On large source-zones, the final stochastic_slip generation step of the above
code might not complete [say if the above script is killed due to exceeding
time limits]. In that case, you can run the unfinished stochastic slip runs
with
[run_unfinished_stochastic_earthquake_tsunami.R](run_unfinished_stochastic_earthquake_tsunami.R).
This will make RDS files with the unfinished results, but does not write them
to the netcdf file (to avoid the chance of multiple files trying to write at
once). Once all the runs have finished, run
[run_merge_additional_stochastic_tsunami_into_netcdf.PBS](run_merge_additional_stochastic_tsunami_into_netcdf.PBS)
to put the results into the netcdf. Finally, you might want to delete the folder R_images_tmp that 
was used to store temporary files. 


