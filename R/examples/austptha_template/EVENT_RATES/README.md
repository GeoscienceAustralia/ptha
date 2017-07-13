This folder contains files to integrate over the tsunami events and nominal
rates which were derived in sub-folders of ../SOURCE_ZONES/

## The main script is:

./compute_station_hazard_curves.R

    Used to compute 'hazard curves' for every station in the model. These are
    written out to a netcdf file [one for each source-zone]

The above script can be run with the PBS script run_compute_rates_and_station_hazard_curves.PBS


## Supplementary codes which are used in the above process are:

./compute_rates_all_sources.R

    used to compute rates for all tsunami events

./make_spatially_variable_source_zone_convergence_rates.R 

    used to map Bird (2003) plate-boundary convergence rates onto the top-edge
    of the unit-sources. 

## Plotting routines

./plot_hazard_curves.R
    Function to plot the wave height exceedance rates for both stochastic slip
    and uniform slip, at a station.
