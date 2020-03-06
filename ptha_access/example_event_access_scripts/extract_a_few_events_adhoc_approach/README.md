
To run the code, you first need to open the file
[extract_time_series.R](extract_time_series.R),
and edit the input variables and the location of the
[get_PTHA_results.R](../../get_PTHA_results.R) script on your machine.

Then, supposing `rptha` and all other add-on packages are installed, and you have
a good internet connection, you should be able to run the code with: 

    Rscript extract_time_series.R

You can also make a quick plot of the results (again, make sure the path to [get_PTHA_results.R](../../get_PTHA_results.R) is correct).

    Rscript quick_plot_stage_at_reference_point.R
