# Examples of scenario data extraction

The scripts in this directory give a few examples of extracting a tsunami
scenarios from the PTHA database.

## Plot tsunami deformation and gauge time-series for a number of events

The code in [gauge_and_deformation_plots](gauge_and_deformation_plots) can be
useful to eyeball many scenarios, to assist with manual event selection.

## Extract a few time-series at some target points, with prescribed maximum-stage values or exceedance-rates, using an ad-hoc approach to get just a few scenarios

See [extract_a_few_events_adhoc_approach](extract_a_few_events_adhoc_approach)

## Initial conditions for scenarios that are similar to historic tsunamis

The code
[./scenarios_similar_to_historical/get_scenarios_similar_to_historical_events.R](./scenarios_similar_to_historical/get_scenarios_similar_to_historical_events.R)
shows how to download and plot the initial conditions for PTHA18 scenarios that
had top-3 goodness of fit when compared with historical data at DART buoys.

There is no guarentee that these scenarios will match the tsunami as observed
elsewhere -- sometimes they will, but not always. In principle one should
expect to get better accuracy using an inverted source that considers data
at your specific site of interest.
