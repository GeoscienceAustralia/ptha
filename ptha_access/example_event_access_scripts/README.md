# Examples of scenario data extraction

The scripts in this directory give a few examples of extracting a tsunami
scenarios from the PTHA database.

## Scenario selection based on exceedance-rates, informed by multiple gauges

See [./multi_site_scenario_selection](./multi_site_scenario_selection)

## Monte-Carlo sampling of scenarios

See [./random_scenarios_non_uniform_and_importance_sampling](./random_scenarios_non_uniform_and_importance_sampling) 
for tutorials on scenario sampling, including the use of stratified and
stratified/importance-sampling, and computation of Monte-Carlo confidence
intervals.

An earlier variant of this tutorial is in
[./random_scenario_sampling](./random_scenario_sampling).

## Plot tsunami deformation and gauge time-series for a number of events

The code in [gauge_and_deformation_plots](gauge_and_deformation_plots) can be
useful to eyeball many scenarios, to assist with manual event selection.

## Extract a few time-series at some target points, with prescribed maximum-stage values or exceedance-rates, using an ad-hoc approach to get just a few scenarios

See [extract_a_few_events_adhoc_approach](extract_a_few_events_adhoc_approach)

## Initial conditions for scenarios that are similar to historic tsunamis

See [./scenarios_similar_to_historical/](./scenarios_similar_to_historical/)

