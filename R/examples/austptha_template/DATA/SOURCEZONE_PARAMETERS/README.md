Key configuration parameters for each source-zone
================================================

The csv file [sourcezone_parameters.csv](sourcezone_parameters.csv) contains
key geometric parameters for each sourcezone, that must be configured correctly
for everything to work. The parameters are:

sourcename
----------

The name of the source. For each 'sourcename', there needs to be a shapefile
named 'sourcename.shp' inside ../SOURCEZONE_CONTOURS, and a file named
'sourcename_downdip.shp' in ../SOURCEZONE_DOWNDIP_LINES. The former file
defines the source's interface geometry (as depth contours).  The latter file defines
the desired location of breaks between unit-sources in the along-strike direction.


segment_name
------------

This must be empty for source-zones without segmentation. Otherwise it gives
the name of the segment, prefixed by an underscore (\_). All segments on the
source-zone must have unique names.

Segmentation allows the use of different source parameters (such as Mw_max, b,
coupling) in different parts of the source-zone. These will affect the rate of
earthquake events that are fully or partly inside the segment. Earthquakes are
always allowed to cross segment boundaries, and earthquakes inside two or more
segments will have their rate determined as a weighted mean of their rates on
each individual segment (the weights are based on the fraction of the
earthquake's unit-sources on each segment).

All source-zones MUST have one (and only one) unsegmented row in the
sourcezone_parameters.csv file. The source-zone may also have a number of
segments, which collectively should cover the entire source, and not overlap
with each other. The relative influence of the unsegmented and segmented
representations can be controlled with the row_weight, mentioned further below. 

This must be finalised before creating the unit-sources using scripts in
[../../SOURCE_ZONES/TEMPLATE/EQ_SOURCE/](../../SOURCE_ZONES/TEMPLATE/EQ_SOURCE/).

approx_unit_source_length
-------------------------

The desired length (along-strike) of unit-sources in km. *Note this is currently
ignored, since the along-strike spacing of unit-sources is determined by the
sourcename_downdip.shp shapefile*. 

This must be finalised before creating the unit-sources using scripts in
[../../SOURCE_ZONES/TEMPLATE/EQ_SOURCE/](../../SOURCE_ZONES/TEMPLATE/EQ_SOURCE/).

approx_unit_source_width
------------------------

The desired width (down-dip) of the unit-sources in km. When discretizing the
source into unit sources, the code will try to make them all this width on average,
although this is constrained by the geometry of the source.  

Note that the width is measured down-dip, not along the surface of the earth.

Note that there should only be one unit_source_width for each source-zone.

This must be finalised before creating the unit-sources using scripts in
[../../SOURCE_ZONES/TEMPLATE/EQ_SOURCE/](../../SOURCE_ZONES/TEMPLATE/EQ_SOURCE/).

rake
----

The desired rake (in degrees) of events on the source. Currently only pure thrust (=90) and pure normal (= -90) are supported.

Note that there should only be one rake for each source-zone.

This must be finalised before creating the unit-sources using scripts in
[../../SOURCE_ZONES/TEMPLATE/EQ_SOURCE/](../../SOURCE_ZONES/TEMPLATE/EQ_SOURCE/).

scaling_relation
----------------

Name of the scaling relation used to relate earthquake area, length and width
to magnitude. It should correspond to a scaling relation accepted by the
`rptha` function `Mw_2_rupture_size`. 

If segmentation is used, the scaling_relation of all segments must be identical
to the scaling_relation on the unsegmented source.

This must be finalised before creating the tsunami events using the scripts in
[../../SOURCE_ZONES/TEMPLATE/TSUNAMI_EVENTS/](../../SOURCE_ZONES/TEMPLATE/TSUNAMI_EVENTS/).

shear_modulus
-------------

Value for the shear modulus (units of 10^10 Pascals). For instance, a value of 3
corresponds to 3x10^10 Pascals. Typical values used on subduction zones range
from 2 to 7, although lower values may occur in shallow areas. 

Note only a single value can be used on each source-zone (i.e. all segments
must have the same value as on the unsegmented source).

This must be finalised before creating the tsunami events using the scripts in
[../../SOURCE_ZONES/TEMPLATE/TSUNAMI_EVENTS/](../../SOURCE_ZONES/TEMPLATE/TSUNAMI_EVENTS/).

segment_boundary_alongstrike_index_lower
----------------------------------------

This should be blank for unsegmented rows. On segments, it gives an integer,
corresponding to the alongstrike-index of the unit-source where the segment
begins. To select this value, you will need to look at the unit_source_grid
shapefile that is created when running
[../../SOURCE_ZONES/TEMPLATE/EQ_SOURCE/produce_unit_sources.R](produce_unit_sources.R). It
will show the locations of the unit sources, and have an attribute table with
their alongstrike and downdip indices.

This must be finalised before computing the event rates in [../../EVENT_RATES/](../../EVENT_RATES).

segment_boundary_alongstrike_index_upper
----------------------------------------

This defines the upper-bound of the alongstrike-index of unit sources that are inside the segment. See
the previous definition for details.

This must be finalised before computing the event rates in [../../EVENT_RATES/](../../EVENT_RATES).

row_weight
----------

This gives the weight of segmented and unsegmented treatments of the
source-zone. All segments should have the same weight. Furthermore, the sum of
the unsegmented weight and the weight assigned to A SINGLE segment should equal one.

For example, if you want to equally weight the unsegmented and segmented
treatments, then you should set the unsegmented weight to 0.5, and all the
segment weights to 0.5. 

If you want to put 99% of the weight on the segmented treatment, then you
should set the unsegmented weight to 0.01, and all the segmented weights to
0.99. 

If you want to put 80% of the weight on the unsegmented treatment, then you
should set the unsegmented weight to 0.8, and all the segmented weights to 0.2.

Row weights cannot be negative! 

This must be finalised before computing the event rates in [../../EVENT_RATES/](../../EVENT_RATES).

use_bird_convergence
--------------------

This takes the value 0 or 1, and should only take a single value on each
source-zone (i.e. segments cannot have a value that differs from the
unsegmented value).

If 0, then use a constant plate convergence rate for the whole source-zone
(provided in the next columns).

If 1, then use the plate boundary convergence data in
'../BIRD_PLATE_BOUNDARIES/sourcezone_traces_table_merged.csv.zip' to define the
relative plate motion rates at each unit source. This is done for each
along-trench unit source (downdip-index = 1) using its nearest point in the
latter file. Note that in this case, depending on options in
[../../EVENT_RATES/config.R](../../EVENT_RATES/config.R), we may either use the
convergent fraction of motion only, or both the convergent fraction and some
part of the transform fraction. This can account for the possibility of
uncertainties in the direction of motion. See the latter config file for more
information.

This must be finalised before computing the event rates in [../../EVENT_RATES/](../../EVENT_RATES).


tectonic_slip
-------------

This should be blank if use_bird_convergence=1. Otherwise, it gives the assumed
plate boundary convergence rate in mm/year, for the entire source. 

This must be finalised before computing the event rates in [../../EVENT_RATES/](../../EVENT_RATES).

convergent_fraction
--------------------

This should be blank if use_bird_convergence=1. Otherwise, it gives the assumed
convergent fraction of the plate boundary convergence rate in mm/year, for the
entire source. 

The implementation is very simple -- the tectonic_slip variable is multiplied
by this before it is used. This means that e.g. (tectonic_slip = 20,
convergent_fraction=0.5) will give identical results to (tectonic_slip = 10,
convergent_fraction=1). 

This must be finalised before computing the event rates in [../../EVENT_RATES/](../../EVENT_RATES).

cmin, cpref, cmax
-----------------

Three values which define the coupling-coefficients for the source-zone or
segment that are used in the logic tree. 

**Beware that these values may or may not be used**, depending on options in
[../../EVENT_RATES/config.R](../../EVENT_RATES/config.R). Parameters in the latter
file also allow for a greater range of coupling values to be used, by interpolating between
the three provided values.

This must be finalised before computing the event rates in [../../EVENT_RATES/](../../EVENT_RATES).

Do not make any of these values equal to zero. If you want zero coupling, then
assign a nonzero value to prob_Mmax_below_Mmin. The reason for this is that we
numerically discretize the coupling coefficients with a logarithmic spacing, so
that in a relative sense we can resolve both small and low couplings. Obviously
one cannot take the log of zero, so zero coupling is not allowed. This should
NOT be interpreted to mean that log(coupling) is uniformly distributed. Rather,
it is a numerical improvement to the discretization. We still use a uniform
distribution to interpolate the prior density among the provided coupling
values - but then we numerically discretize that density with a log spacing.

cmin_p,cpref_p,cmax_p
---------------------

**NOT USED.** Three values which define the weights of coupling-coefficients
for the source-zone or segment that are used in the logic tree. 

*BEWARE CURRENTLY THESE VALUES ARE IGNORED, AND EACH COUPLING COEFFICIENT IS
ASSIGNED THE SAME WEIGHT*.


bmin, bpref, bmax
-----------------

Three values which define the Gutenberg-Richter b values for the source-zone or
segment that are used in the logic tree. 

This must be finalised before computing the event rates in [../../EVENT_RATES/](../../EVENT_RATES).

bmin_p,bpref_p,bmax_p
---------------------

**NOT USED.** Three values which define the weights of Gutenberg-Richter b values
for the source-zone or segment that are used in the logic tree. 

*BEWARE CURRENTLY THESE VALUES ARE IGNORED, AND EACH b COEFFICIENT IS
ASSIGNED THE SAME WEIGHT*. This can be changed by modifying code in
[../../EVENT_RATES/compute_rates_all_sources.R](../../EVENT_RATES/compute_rates_all_sources.R)

mw_max_observed
---------------

The largest earthquake thought to have historically occurred on the source-zone
or segment. This provides a lower-limit to logic-tree values of Mw_max with
non-zero weight. Note that a small number is added to this before use to ensure
that the largest historical earthquake is possible under our logic-tree
branches. See
[../../EVENT_RATES/config.R](../../EVENT_RATES/config.R)

This must be finalised before computing the event rates in [../../EVENT_RATES/](../../EVENT_RATES).


prob_Mmax_below_Mmin
----------------------

The prior probability that the source zone does NOT have any seismicity above our Mw-min.
This can be used to treat source-zones that might be inactive. 

notes
------

Comments, which are ignored for computation purposes
