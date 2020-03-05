# Examples of scenario data extraction

The scripts in this directory give a few examples of extracting a tsunami
scenarios from the PTHA database.

## Extracting a few time-series at some target points, with prescribed maximum-stage values or exceedance-rates, using an ad-hoc approach to get just a few scenarios

To run the code, you first need to open the file
[./extract_a_few_events_adhoc_approach/extract_time_series.R](./extract_a_few_events_adhoc_approach/extract_time_series.R),
and edit the input variables and the location of the
[get_PTHA_results.R](../get_PTHA_results.R) script on your machine.

Then, supposing `rptha` and all other add-on packages are installed, and you have
a good internet connection, you should be able to run the code with: 

    Rscript extract_time_series.R

You can also make a quick plot of the results (again, make sure the path to [get_PTHA_results.R](../get_PTHA_results.R) is correct).
    Rscript quick_plot_stage_at_reference_point.R

## Initial conditions for scenarios that are similar to historic tsunamis

The code
[./scenarios_similar_to_historical/get_scenarios_similar_to_historical_events.R](./scenarios_similar_to_historical/get_scenarios_similar_to_historical_events.R)
shows how to download and plot the initial conditions for PTHA18 scenarios that
had top-3 goodness of fit when compared with historical data at DART buoys.

There is no guarentee that these scenarios will match the tsunami as observed
elsewhere -- sometimes they will, but not always. In principle one should
expect to get better accuracy using an inverted source that considers data
at your specific site of interest.
