This folder contains files to integrate over the tsunami events and nominal
rates which were derived in sub-folders of ../SOURCE_ZONES/. 

The latter codes, which relate to making tsunami events on source zones, must
have been run before trying to run codes here.


## Running 

To run this code, you must have created the unit-sources, run the tsunami models,
and created the earthquake events. See sub-folders in
[../SOURCE_ZONES/TEMPLATE/](../SOURCE_ZONES/TEMPLATE) for more information.

You also need to decide on 'final' parameter values for the source-zone, and how
it is segmented. These are specified in
[../DATA/SOURCEZONE_PARAMETERS/sourcezone_parameters.csv](../DATA/SOURCEZONE_PARAMETERS/sourcezone_parameters.csv).
For more information on these parameters, see 
[../DATA/SOURCEZONE_PARAMETERS/sourcezone_parameters/README.md](../DATA/SOURCEZONE_PARAMETERS/README.md).

You also need to edit [config.R](config.R) in this folder, to set parameters
used in the event rate computation. See the comments in that script for more details.
*It is essential that users understand what these parameters do, otherwise your computation
might be completely different to what you wanted!*

One input to the [config.R](config.R) script consists of a csv file describing
functions used bias-correct event weights (with the same magnitude and general
location) based on their peak slip. We used the script
[event_properties_and_GOF.R](event_properties_and_GOF.R) to make these, based
on the DART buoy comparisons. 

Then, assuming everything is correctly configured, the codes can be run (on
NCI) with:

    qsub run_compute_rates_and_station_hazard_curves.PBS

With slight modifications, one could also just 'source' the above script on
another machine. 

Once that is done, you should also do the following:

1. Run [create_nc_file_peak_stage_unit_sources.R](create_nc_file_peak_stage_unit_sources.R) to make fast-access information for plotting. 
2. To make easy-access output files you need to run [tsunami_stage_exceedance_rates_to_csv.R](tsunami_stage_exceedance_rates_to_csv.R). 
3. To make a shapefile for 'general-plot' purposes, run [clean_shapefiles_for_plotting.R](clean_shapefiles_for_plotting.R).

You can also use [quick_station_stage_exceedance_rates.R](quick_station_stage_exceedance_rates.R) to make deaggregation plots on-demand.

## Details

The main scripts are [compute_station_hazard_curves.R](compute_station_hazard_curves.R)
and [compute_rates_all_sources.R](compute_rates_all_sources.R)

[compute_rates_all_sources.R](compute_rates_all_sources.R) is used to compute
rates for all tsunami events on all source-zones, using information on 
plate tectonic motions and historical earthquakes. 

[compute_station_hazard_curves.R](compute_station_hazard_curves.R) combines the
above information with the tsunami propagation results to compute 'hazard
curves' for every station in the model. These are written out to a netcdf file
[one for each source-zone]. *If you have many and/or large source-zones, this
code may take a long time to run (e.g. around 36 hours for the final PTHA18
runs, on a 16-core node of raijin). To work around this, you can optionally run
the code on a subset of source-zones only, by passing an integer commandline
argument ranging from 1 to 6 (details are currently hard-coded for the PTHA18
source-zones).  This means the 'batches' can be run consecutively on separate
nodes, which may speed up the overall run progress. The general idea is that
for each source-zone we call a function* `check_if_file_is_ok` *to see whether
that source-zone should be run.  This can be manipulated to only run a subset
of source-zones. See the source code for further information.*

Supplementary codes which are used in the above process are:

[append_variable_mu_variables_to_event_netcdf.R](append_variable_mu_variables_to_event_netcdf.R), used
to insert new columns into the netcdf files to store some results under the
assumption of variable shear modulus (this is only different to the
constant-shear-modulus case for thrust events). In principle it would
be nicer to add these variables when the files are originally created. But
in reality, this functionality was 'bolted on' to the PTHA code after 
we were already dealing with the constant shear modulus case, which shows up
in some of the design

[make_spatially_variable_source_zone_convergence_rates.R](make_spatially_variable_source_zone_convergence_rates.R),
used to map Bird (2003) plate-boundary convergence rates onto the top-edge of
the unit-sources. 

[integrated_rate_given_stage.R](integrated_rate_given_stage.R)
used to compute the stage-vs-exceedance rate curves at each hazard point,
integrated over all source zones. Note that for the credible interval
computation, we assume that uncertainties at all source-zones behave
co-monotonically (i.e. there is no cancellation of uncertainty due to
averaging epistemic uncertainties at different sites). The co-mononotonic
assumption follows Davies et al (2017), and is a conservative treatment,
used to avoid having to specify the dependency structure of our epistemic
uncertainties across source-zones. 

[gcmt_subsetter.R](gcmt_subsetter.R) was used to extract subsets of the GCMT data
for later analysis.

[back_calculate_convergence.R](back_calculate_convergence.R) contains code to 
calculate and plot the slip rate on each unit source implied by the model. This
can be compared with the plate-tectonic models which were used to determine the
event rates. In combination with a numerical optimization routine, it may also
be used to determine how much the conditional probability of source-zone-edge
earthquake events should be increased in order to best match the spatial distribution
of plate convergence.

[check_event_netcdf_files.R](check_event_netcdf_files.R) performs various
sanity checks on the netcdf files that are updated by
[compute_rates_all_sources.R](compute_rates_all_sources.R). 

[event_properties_and_GOF.R](event_properties_and_GOF.R) shows some investigation of relationships
between rupture geometry and goodness-of-fit. This is used to create 'bias
adjustments' for heterogeneous and variable-uniform slip events, by comparing the statistical
properties of 'modelled events which are similar to
deep-ocean-tsunami-observations' to 'all modelled events'. 

## Plotting routines

[quick_station_stage_exceedance_rates.R](quick_station_stage_exceedance_rates.R)
contains code to make diagnostic plots that help to understand the result at a single point.
The plots are designed both to convey the results, and try to highlight any
'problematic' issues (e.g. related to convergence of the hazard results / cross-checks of
results between different files, etc). It uses
[plot_hazard_curves_utilities.R](plot_hazard_curves_utilities.R) which is has
lots of useful disgnostic plots that can help to understand the analysis. 

## Other routines

[tsunami_stage_exceedance_rates_to_csv.R](tsunami_stage_exceedance_rates_to_csv.R)
Script to convert the netcdf files containing stage-vs-return-period curves,
into pointwise values of stages at a range of return periods, in a csv and
shapefile format. The resulting outputs have some rounding, and contain less
information than the netcdf output. However they are provided because many users
find these formats easy to work with.

[stage_range_summary.R](stage_range_summary.R) contains basic checks on the performance
of the uniform, stochastic, and variable uniform models vs DART measurements. 

[slip_simulator.R](slip_simulator.R) is a simple (and experimental) code to
place fixed-size earthquake events on a source-zone, so as to spatially match
some idea of moment conservation.

[working_with_rate_curves.R](working_with_rate_curves.R) gives an example of working
with the rate curves for a particular source-zone.

[variable_mu_checks.R](variable_mu_checks.R) makes some plots to compare Mw-frequency curves
with fixed and variable mu. NOTE: This code assumes that the netcdf files have
NOT already had bias adjustments applied (i.e. if you run them on the final
output files, you'll get strange results!). The were used to help explore
various bias adjustment methods during preliminary stages of the analysis.

[convergence_rates_plots_comparison_methods.R](convergence_rates_plots_comparison_methods.R) is
now outdated, but was used to highlight the significance of edge-correction of
event rates. FIXME: It needs to be updated to include the peak slip limitation
that was not accounted for in earlier runs.

[earthquake_rate_comparisons.R](earthquake_rate_comparisons.R) is used to compute the modelled
scenario rates in various regions.

*(DEFUNCT)* [plot_peak_stage_1m_slip.R](plot_peak_stage_1m_slip.R) is an old
plotting code, which is superceeded by
[quick_station_stage_exceedance_rates.R](quick_station_stage_exceedance_rates.R)
which makes plots of similar information.

*(DEFUNCT)* [plot_hazard_curves.R](plot_hazard_curves.R) is an old plotting
routine, which is superceeded by
[quick_station_stage_exceedance_rates.R](quick_station_stage_exceedance_rates.R) 

