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

Template codes used for a 'real' large-scale PTHA for Australia. Includes code
for generating unit-source tsunami, running hydrodynamic models, creating
stochastic and uniform slip events, and hazard computation.
