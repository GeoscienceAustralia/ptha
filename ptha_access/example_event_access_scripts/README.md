# Examples of scenario data extraction

The code in this directory gives examples of extracting a few tsunami events
from the PTHA database.

## Extracting time-series at some target points, with given rates or stages.

To run the code, you first need to open the file 'extract_time_series.R', and edit
the path of the script which it sources, named 'get_PTHA_results.R'. 

You should ensure it points to the location of that script on your machine.

Then, supposing `rptha` and all other add-on packages are installed, and you have
a good internet connection, you should be able to run the code with: 

    Rscript extract_time_series.R
    Rscript quick_plot_stage_at_reference_point.R

## Initial conditions for scenarios that are similar to historic tsunamis

The code get_scenarios_similar_to_historical_events.R shows how to download and
plot the initial conditions for PTHA18 scenarios that had top-3 goodness of fit
when compared with historical data at DART buoys.
