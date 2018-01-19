This code gives an example of extracting a few tsunami events from the PTHA database

To run the code, you first need to open the file 'extract_time_series.R', and edit
the path of the script which it sources, named 'get_PTHA_results.R'. 

You should ensure it points to the location of that script on your machine.

Then, supposing `rptha` and all other add-on packages are installed, and you have
a good internet connection, you should be able to run the code with: 

    Rscript extract_time_series.R
    Rscript quick_plot_stage_at_reference_point.R
