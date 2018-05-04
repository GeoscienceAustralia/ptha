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

Then, assuming everything is correctly configured, the codes can be run (on
NCI) with:

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

## Plotting routines

[plot_hazard_curves.R](plot_hazard_curves.R) contains a function to plot the
wave height exceedance rates for both stochastic slip and uniform slip, at a
station. There are lots of useful disgnostic plots here that can help to understand the analysis.

[quick_station_stage_exceedance_rates.R](quick_station_stage_exceedance_rates.R)
contains code to 'independently' calculate the stage exceedance rate curves at a point
(i.e. not relying on the above codes [compute_station_hazard_curves.R](compute_station_hazard_curves.R)
or  [integrated_rate_given_stage.R](integrated_rate_given_stage.R) ). This was done
as a partial check on the correctness of the latter routines (in the
absence of more exact tests, we can at least confirm that we get the same answer
using 2 different methods -- that's the idea behind this routine).

[plot_peak_stage_1m_slip.R](plot_peak_stage_1m_slip.R) can be used to plot the peak stage
at a chosen hazard point resulting from every single unit-source tsunami. It's
useful as a sanity check on the tsunami model results.

[convergence_rates_plots_comparison_methods.R](convergence_rates_plots_comparison_methods.R) can
be used to graphically compare back-calculated convergence given different assumptions. The main
purpose is to highlight the significance of edge-correction of event rates.

## Other routines

[stage_range_summary.R](stage_range_summary.R) contains basic checks on the performance
of the uniform, stochastic, and variable uniform models vs DART measurements. 

[slip_simulator.R](slip_simulator.R) is a simple (and experimental) code to
place fixed-size earthquake events on a source-zone, so as to spatially match
some idea of moment conservation.

