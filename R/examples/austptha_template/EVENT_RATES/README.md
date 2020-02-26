# Codes for the scenario rate computation in PTHA18

This folder contains scripts to do the source-zone rate computations in PTHA18,
distribute the rates over scenarios, compute the hazard curves, etc. 

## Running 

To run the scenario rate calculations, you must have created the unit-sources,
run the tsunami models, and created the earthquake events. See sub-folders in
[../SOURCE_ZONES/TEMPLATE/](../SOURCE_ZONES/TEMPLATE) for more information.

You also need to decide on 'final' parameter values for the source-zone, and how
it is segmented. These are specified in
[../DATA/SOURCEZONE_PARAMETERS/sourcezone_parameters.csv](../DATA/SOURCEZONE_PARAMETERS/sourcezone_parameters.csv).
For more information on these parameters, see 
[../DATA/SOURCEZONE_PARAMETERS/sourcezone_parameters/README.md](../DATA/SOURCEZONE_PARAMETERS/README.md).

You also need to edit [config.R](config.R) in this folder, to set parameters
used in the event rate computation. See the comments in that script for more
details.

One input to the [config.R](config.R) script consists of a csv file describing
functions used bias-correct event weights (with the same magnitude and general
location) based on their peak slip. We used the script
[event_properties_and_GOF.R](event_properties_and_GOF.R) to make these, based
on the DART buoy comparisons. 

Then, assuming everything is correctly configured, the codes can be run (on
NCI) with:

    qsub run_compute_rates_and_station_hazard_curves.PBS

With slight modifications, one could also just 'source' the above script on
another machine. Alternatively, you could look at the sequence of commands in
the script and then run them some other way (e.g. this can allow the rate curve
calculations to be split over multiple nodes).

Once that is done, you should also do the following to update various derived
results (*NOTE: you might need to edit output paths in these scripts first, they
are not tightly integrated with the previous computations*):

1. Run [create_nc_file_peak_stage_unit_sources.R](create_nc_file_peak_stage_unit_sources.R) to make fast-access information for plotting. Note this is unaffected by the scenario rates, so it only needs to be run when the tsunami-unit-sources are updated.
2. To make easy-access output files you need to run [tsunami_stage_exceedance_rates_to_csv.R](tsunami_stage_exceedance_rates_to_csv.R). 
3. To make a shapefile for 'general-plot' purposes, run [clean_shapefiles_for_plotting.R](clean_shapefiles_for_plotting.R).

You can also use
[revised_quick_station_stage_exceedance_rates.R](revised_quick_station_stage_exceedance_rates.R)
to make 'standard station deaggregation plots' on-demand. To assist with batch-production of such
plots, see [revised_quick_station_plots_all_sites.R](revised_quick_station_plots_all_sites.R)
which runs a group of hazard points based on their longitude.

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
source-zones - see the function `skip_this_source_zone` for details). This means
multiple subsets can be submitted at once on separate nodes, which may speed up
the overall run progress.*

## Supplementary codes which are used in the above process are:

[config.R](config.R) controls key run variables.

[append_variable_mu_variables_to_event_netcdf.R](append_variable_mu_variables_to_event_netcdf.R)
is used to insert new columns into the netcdf files. The new columns store some
event rates under the assumption of variable shear modulus (for thrust events).
In principle it would be nicer to add these variables when the files are
originally created. But in reality, this functionality was 'bolted on' to the
PTHA code after we were already dealing with the constant shear modulus case.
Hence the need for this script.

[make_spatially_variable_source_zone_convergence_rates.R](make_spatially_variable_source_zone_convergence_rates.R),
used to map Bird (2003) plate-boundary convergence rates onto the top-edge of
the unit-sources. 

[integrated_rate_given_stage.R](integrated_rate_given_stage.R) is
used to compute the stage-vs-exceedance rate curves at each hazard point,
integrated over all source zones. Note that for the credible interval
computation, we assume that uncertainties at all source-zones behave
co-monotonically (i.e. there is no cancellation of uncertainty due to
averaging epistemic uncertainties at different sites). The co-mononotonic
assumption follows Davies et al (2017), and is a conservative treatment,
used to avoid having to specify the dependency structure of our epistemic
uncertainties across source-zones. 

[gcmt_subsetter.R](gcmt_subsetter.R) is used to extract subsets of the GCMT data
on source-zones. 

[back_calculate_convergence.R](back_calculate_convergence.R) contains code to 
calculate and plot the slip rate on each unit source implied by the model. This
can be compared with the plate-tectonic models which were used to determine the
event rates. In combination with a numerical optimization routine, it may also
be used to determine how much the conditional probability of source-zone-edge
earthquake events should be increased in order to best match the spatial
distribution of plate convergence.

[check_event_netcdf_files.R](check_event_netcdf_files.R) performs various
sanity checks on the netcdf files that are updated by
[compute_rates_all_sources.R](compute_rates_all_sources.R). 

[event_properties_and_GOF.R](event_properties_and_GOF.R) shows some
investigation of relationships between rupture geometry and goodness-of-fit.
This was used to create 'bias adjustments' for heterogeneous and
variable-uniform slip events, by comparing the statistical properties of
'modelled events which are similar to deep-ocean-tsunami-observations' to 'all
modelled events'. It creates files which are referenced in [config.R](config.R).
Note you have to manually edit the script to specify whethere variable shear modulus is used.


## Plotting routines are

*(NOTE -- it is preferable to use the revised versions of the routines in this paragraph-- see the section below named "Updated stage-vs-exceedance-rate percentile uncertainty calculations"))* 
[quick_station_stage_exceedance_rates.R](quick_station_stage_exceedance_rates.R)
contains code to make diagnostic plots that help to understand the result at a
single point.  The plots are designed both to convey the results, and try to
highlight any 'problematic' issues (e.g. related to convergence of the hazard
results / cross-checks of results between different files, etc). It uses
[plot_hazard_curves_utilities.R](plot_hazard_curves_utilities.R) which is has
lots of useful disgnostic plots that can help to understand the analysis. To
run many points in one go, see the code
[quick_station_plots_all_sites.R](quick_station_plots_all_sites.R), which can
run all points (or just a set of points with nearby longitudes). To do the
latter on a a shared-memory parallel machine, see
[run_quick_station_plots_on_subsets_of_sites.PBS](run_quick_station_plots_on_subsets_of_sites.PBS).

[event_dart_coverage_vs_distance.R](event_dart_coverage_vs_distance.R) can make
a plot showing the observed stage-ranges as a percentile of the corresponding
family of model scenarios. Before this is run, you must haev run
[event_properties_and_GOF.R](event_properties_and_GOF.R).

## Other routines

[tsunami_stage_exceedance_rates_to_csv.R](tsunami_stage_exceedance_rates_to_csv.R)
Script to convert the netcdf files containing stage-vs-return-period curves,
into pointwise values of stages at a range of return periods, in a csv and
shapefile format. The resulting outputs have some rounding, and contain less
information than the netcdf output. However they are provided because many users
find these formats easy to work with.

[stage_range_summary.R](stage_range_summary.R) contains checks on the
performance of the uniform, stochastic, and variable uniform models vs DART
measurements (including hypothesis testing relating to the median-coverage-statistic)

[earthquake_rate_comparisons_PTHA18_report.R](earthquake_rate_comparisons_PTHA18_report.R) was used to
compute the modelled scenario rates in various regions in the 
[PTHA18 report](http://dx.doi.org/10.11636/Record.2018.041) (for comparison with
regional estimates in the literature)

[earthquake_rate_comparisons_PAGEOPH_paper.R](earthquake_rate_comparisons_PAGEOPH_paper.R) was used to
compute the modelled scenario rates in various regions in [this
paper](https://link.springer.com/article/10.1007%2Fs00024-019-02299-w). As compared
with the aforementioned version, there are some changes to the
region-for-comparison and the interpolation method.

## Exploratory routines (not directly used elsewhere)

[slip_simulator.R](slip_simulator.R) is a simple (and experimental) code to
greedily place fixed-size earthquake events on a source-zone, so as to spatially match
some idea of moment conservation. This is not used directly in the above analysis,
but it independently shows that edge-effects should be accounted for.

[variable_mu_checks.R](variable_mu_checks.R) makes some plots to compare
Mw-frequency curves with fixed and variable mu. NOTE: This code assumes that
the netcdf files have NOT already had bias adjustments applied (i.e. if you run
them on the final output files, you'll get strange results!). The were used to
help explore various bias adjustment methods during preliminary stages of the
analysis.

[convergence_rates_plots_comparison_methods.R](convergence_rates_plots_comparison_methods.R) 
shows long-term convergence rates implied different conditional probability approaches.

*(DEFUNCT)* [plot_peak_stage_1m_slip.R](plot_peak_stage_1m_slip.R) is an old
plotting code, which is superceeded by
[quick_station_stage_exceedance_rates.R](quick_station_stage_exceedance_rates.R)
which makes plots of similar information.

*(DEFUNCT)* [plot_hazard_curves.R](plot_hazard_curves.R) is an old plotting
routine, which is superceeded by
[quick_station_stage_exceedance_rates.R](quick_station_stage_exceedance_rates.R) 

[working_with_rate_curves.R](working_with_rate_curves.R) gives an example of working
with the rate curves for a particular source-zone.

## Updated stage-vs-exceedance-rate percentile uncertainty calculations

These scripts apply some revised stage-vs-exceedance-rate percentile
uncertainty calculations, as described in [Section 3.5 of this paper](https://link.springer.com/article/10.1007/s00024-019-02299-w).

The percentiles of the exceedance-rates are thus updated compared with the
original PTHA18 report. The equations for the hazard calculation imply that the
logic-tree mean exceedance-rate is unaffected. However our revised
implementation induced minor numerical changes even in the mean exceedance-rate
(typically in high significant figures -- not large enough to show up in
plots). This was due to a change in interpolation at one point in the
calculation. It is just as valid as the old approach, but much more
computationally efficient in the context of the revised approach. 

[revised_station_hazard_curves_PREPROCESSING.R](revised_station_hazard_curves_PREPROCESSING.R) preprocesses source-zone logic-tree branches to support the calculation.

[revised_station_hazard_curves.R](revised_station_hazard_curves.R) computes stage-vs-exceedance-rate percentile curves. This version of the script uses logic-tree sampling to reduce the computational effort, and only works for stochastic-slip + fixed rigidity. For the full calculation, instead see the script below.

[revised_station_hazard_curves_FINAL.R](revised_station_hazard_curves_FINAL.R) computes stage-vs-exceedance-rate percentile curves. This is very similar to [revised_station_hazard_curves.R](revised_station_hazard_curves.R), but it doesn't use logic-tree sampling to reduce the computational effort, and applies to all combinations of slip/rigidity model. So it is much more computationally demanding. For this reason, the script takes input arguments that allow specifying some subset of points. That enables a distributed-computing approach to the calculation.

[submit_all_PBS_revised_station_hazard_curves_FINAL.R](submit_all_PBS_revised_station_hazard_curves_FINAL.R) creates 50 PBS scripts which do the calculations for [revised_station_hazard_curves_FINAL.R](revised_station_hazard_curves_FINAL.R) on 50 nodes (by splitting the points up). This enabled the calculation to be done in a reasonable length of time.

[revised_station_hazard_curves_FINAL_MERGE.R](revised_station_hazard_curves_FINAL_MERGE.R) combines the results created using [submit_all_PBS_revised_station_hazard_curves_FINAL.R](submit_all_PBS_revised_station_hazard_curves_FINAL.R), producing files that look the same as if they were created on a single node.

[revised_ari500_station_hazard_curves_extract.R](revised_ari500_station_hazard_curves_extract.R) processes outputs from the previous script to compute the maximum-stage at ARI=500.

[revised_station_hazard_curves_MAKE_NETCDF.R](revised_station_hazard_curves_MAKE_NETCDF.R) makes netcdf files with the stage-vs-exceedance-rate curves, based on these revised calculations. The file format is the same as that which comes from [compute_station_hazard_curves.R](compute_station_hazard_curves.R), but the filename has `revised1_` prepended.

[revised_integrated_rate_given_stage.R](revised_integrated_rate_given_stage.R) make a netcdf file with the integrated stage-vs-exceedance-rate calculations. It is similar to [integrated_rate_given_stage.R](integrated_rate_given_stage.R), but uses the revised calculations, and the output filename has `revised1_` prepended.

[revised_tsunami_stage_exceedance_rates_to_csv.R](revised_tsunami_stage_exceedance_rates_to_csv.R) make some derivative products, similar to [tsunami_stage_exceedance_rates_to_csv.R](tsunami_stage_exceedance_rates_to_csv.R), but uses the revised calculations, with filenames having `revised1_` prepended.

[revised_clean_shapefiles_for_plotting.R](revised_clean_shapefiles_for_plotting.R) produces some shapefiles for convenient plotting. It is just like [clean_shapefiles_for_plotting.R](clean_shapefiles_for_plotting.R) but using the revised stage percentile results.

[revised_quick_station_stage_exceedance_rates.R](revised_quick_station_stage_exceedance_rates.R) updates the earlier plotting code to use the newer percentiles (i.e. [quick_station_stage_exceedance_rates.R](quick_station_stage_exceedance_rates.R)).

[revised_quick_station_plots_all_sites.R](revised_quick_station_plots_all_sites.R) updates the earlier plotting code to use the newer percentiles (i.e. [quick_station_plots_all_sites.R](quick_station_plots_all_sites.R)).
