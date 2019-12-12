Examples
========

This folder contains various examples of using the rptha package. The
scripts can be used to see how the package works, and may be adapted to
new applications.

source_contours_2_unit_sources
------------------------------

This contains
[example_code](source_contours_2_unit_sources/produce_unit_sources.R) to make
tsunami unit sources from source contours and a
[tutorial](source_contours_2_unit_sources/tutorial.md) on its usage.

event_rates
-----------

Contains an [example script](event_rates/single_source_rate_computation.R) to
compute the rate of earthquake events on a source.


combine_tsunami_sources
-----------------------

[Example script](combine_tsunami_sources/combine_tsunami_sources.R) to make
random slip tsunami from unit source initial conditions.


make_hazard_points
------------------

Although not heavily using the rptha package, this contains example scripts to
[make hazard points](make_hazard_points/make_hazard_pts.R) (i.e. offshore
points where the tsunami tide-gauges are recorded)


austptha_template
-----------------

Template codes used for the 2018 Australian Probabilistic Tsunami Hazard
Assessment (PTHA18). Includes code for generating unit-source tsunami, running
hydrodynamic models, creating stochastic and uniform slip events, and hazard
computation.


okada_displacements_ptha18_scenarios
-------------------------------------

Code which uses the 2018 Australian PTHA geometries to compute the unit-source
easting/northing/vertical Okada displacements. This is modified from code used
in the Australian PTHA
([here](./austptha_template/SOURCE_ZONES/TEMPLATE/EQ_SOURCE/). Modification was
necessary because the latter only stored the Kajiura-smoothed vertical
displacement.
