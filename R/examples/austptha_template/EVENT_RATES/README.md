This folder contains files to integrate over the tsunami events and nominal
rates which were derived in sub-folders of ../SOURCE_ZONES/. 

The latter codes, which relate to making tsunami events on source zones, must
have been run before trying to run codes here.


## Running 

Once the prerequesite steps are complete, edit [config.R](config.R) in this
folder to set parameters used in the event rate computation. Then, assuming
everything is correctly configured, the codes can be run (on NCI) with:

    qsub run_compute_rates_and_station_hazard_curves.PBS

With slight modifications, one could also just 'source' the above script on
another machine.



## Details

The main script is [compute_station_hazard_curves.R](compute_station_hazard_curves.R)

It is used to compute 'hazard curves' for every station in the model. These are
written out to a netcdf file [one for each source-zone]

Supplementary codes which are used in the above process are:

[compute_rates_all_sources.R](compute_rates_all_sources.R), used to compute
rates for all tsunami events

[make_spatially_variable_source_zone_convergence_rates.R](make_spatially_variable_source_zone_convergence_rates.R),
used to map Bird (2003) plate-boundary convergence rates onto the top-edge of
the unit-sources. 

## Plotting routines

[plot_hazard_curves.R](plot_hazard_curves.R) contains a function to plot the
wave height exceedance rates for both stochastic slip and uniform slip, at a
station.
