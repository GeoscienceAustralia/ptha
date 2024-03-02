# Scripts to process the original data, make figures, and produce post-processed outputs

**You do not need to run the code to get the outputs (they are provided in the [../post-processed](../post-processed) folder).**

This folder provides the codes that were used to process the data. 

If you just want to use the data, please do not re-run the code. Some of the original data was removed prior to archiving, and attempts to run the codes may inadvertantly delete some post-processed datasets!

## Required software 

The analysis was run on a linux machine (Ubuntu 22.04). The processing scripts are all written in the R language. They were run with R version 4.2.2, but are likely to work with other versions of R. 

The R interpreter can be [downloaded here](https://cran.r-project.org/). R packages can be installed from inside R with the `install.packages` command. To get the packages used here, the command is:
```r
install.packages(c('sp', 'geosphere', 'fields', 'RColorBrewer', 'magick'))
```
The analysis will optionally use the following packages, if available: 
```r
install.packages(c('fftw', 'maptools'))
```
But they are not essential, and may be harder to install.

### TPXO interface

The folder [TPXO9_interface](TPXO9_interface) contains R code to interface with a local installation of the tidal-prediction software TPX09. It needs to be setup specifically for the user's machine. See README therein for information on how to do this.

## Main scripts to implement the analysis

After the TPXO interface was setup, the analysis was run with the following command (in a bash shell):
```
source run_all_calculations.sh
```

This calls the following scripts:

* [plot_pressure_gauge_locations_and_write_metadata_to_csv.R](plot_pressure_gauge_locations_and_write_metadata_to_csv.R) - read the pressure gauge data, plot the locations, and write the locations to a csv.
* [plot_tide_gauge_locations_and_write_metadata_to_csv.R](plot_tide_gauge_locations_and_write_metadata_to_csv.R) - read the tide gauge data, plot the locations, and write the locations to a csv.
* [plot_pressure_time_series_and_isolate_short_period_waves.R](plot_pressure_time_series_and_isolate_short_period_waves.R) - Plot each MSLP gauge in a number of ways for QC purposes, extract waves with period less than 2 hours, and export results to csv files.
* [plot_tide_gauge_time_series_and_isolate_short_period_waves.R](plot_tide_gauge_time_series_and_isolate_short_period_waves.R) - Plot each tide gauge in a number of ways for QC purposes. Manually fix some records. De-tide each record, and export results to csv files.
* [plot_gauges_near_each_other.R](plot_gauges_near_each_other.R) - Plot the results for pairs of gauges that are located close to each other. This can illustrate the effect of different instruments and different down-sampling methods.
* [compute_gauge_temporal_interval.R](compute_gauge_temporal_interval.R) - Computes some statistics about the gauge sampling rate and duration.
* [plots_for_paper.R](plots_for_paper.R) - Makes some additional plots used in the paper.
* [merge_some_figures.R](merge_some_figures.R) - Combines some of the plots into single panels for the paper.

The above scripts make use of these helper routines:

* [global_variables.R](global_variables.R) - Variables that are useful in many scripts.
* [get_simple_world_map_data.R](get_simple_world_map_data.R) - Polygon data for a world map. The structure is designed to work around future deprecation of the R `maptools` package.
* [spectral_highpass_filter.R](spectral_highpass_filter.R) - Function to separate a time-series into high-frequency and low-frequency components.
* [detiding.R](detiding.R) - Utilities for detiding the tide-gauge data.
* [parse_gauge_data.R](parse_gauge_data.R) - Functions to parse the original pressure and tide-gauge data into a consistent format.
* [create_README_in_postprocessed_folder.R](create_README_in_postprocessed_folder.R) - Make documentation in the postprocessed folder (avoiding manual updates when we change things).
* [create_README_in_postprocessed_graphical_checks_folder.R](create_README_in_postprocessed_graphical_checks_folder.R) - Make documentation in the postprocessed graphical checks folder (avoiding manual updates when we change things).

as well as the binary data [wrld_simpl.RDS](wrld_simpl.RDS) and codes in the [TPXO9_interface](TPXO9_interface) folder.


