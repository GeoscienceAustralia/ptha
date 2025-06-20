# Generate random scenarios

The random scenario generation was done in separate batches over time, due to
of limitations in our quarterly supercomputer quota. 

Each batch is in a subfolder like `set_range_of_mw_and_centroi*`. These folders
contain scripts that do the sampling which are almost identical, except for the:
* Historical events that inform the sampling
* Number of scenarios for each historical event and slip model type
  * In combination we always have 60 scenarios per event and slip model, but these were done as 2 sets of 30 scenarios in the first and second batch of runs, and in sets of 60 thereafter. 

## Code in this folder
* [get_parameters_for_events.R](get_parameters_for_events.R) extracts summary information for the historical events. Some scripts below rely on this having been run.
* [get_target_scenario_focal_mechanism.R](get_target_scenario_focal_mechanism.R) extracts focal mechanism info from GCMT or ISC-GEM. Note this is ONLY used for plotting, not specifically for defining the random scenarios that we model.
* [get_alongstrike_range_for_plotting.R](get_alongstrike_range_for_plotting.R) also extracts information that is helpful for plotting.
* [plot_events.R](plot_events.R) makes a plot.

## Download the scenario datasets

This folder should contain the random scenarios as well as the code used to create them.

However only the code is provided here because the resulting files are large. Thus they are available for download separately at: https://thredds.nci.org.au/thredds/fileServer/fj6/PTHA/Nearshore_testing_2025/ptha18_scenarios_random.tar.bz2

The above tar.bz2 file can be extracted (e.g. `tar -jxf ptha18_scenarios_random.tar.bz2`) to produce a folder named `ptha18_scenarios_random` (i.e. the same name as this folder) which contains sub directories with the same names as subdirectories here. They can merged with the current folder. 
