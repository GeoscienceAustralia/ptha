# **Guide to accessing the 2018 Australian Probabilistic Tsunami Hazard Assessment (PTHA18) results**

# ***NOTE: This is a beta release. Results may be updated if any issues are identified. Please report any problems via the github issues page, or send an email to the maintainer or to hazards@ga.gov.au***

This guide explains how to access basic tsunami hazard information in easy-to-use csv and
shapefile formats. 

For access to detailed PTHA18 information, you should instead consult the [detailed README](DETAILED_README.md), which explains access to:
* scenario earthquakes
* scenario tsunami initial conditions
* scenario wave time-series at 'hazard points'
* additional average-return-interval calculations

The methodologies used in this study, and associated testing, are discussed in
[the project report](http://dx.doi.org/10.11636/Record.2018.041) and 
[this GJI publication](https://doi.org/10.1093/gji/ggz260) and 
[this PAGEOPH publication](https://link.springer.com/article/10.1007/s00024-019-02299-w)
and [this talk](https://www.youtube.com/watch?v=brRy6YjwnlA&list=PL0jP_ahe-BFk3499UEvm-YdlOiE9HFbCg&index=4&t=0s).
Code used to conduct the analysis is available open-source
in the [ptha package](https://github.com/GeoscienceAustralia/ptha). This includes
generic software for PTHA, and also a folder with [project specific
scripts](https://github.com/GeoscienceAustralia/ptha/tree/master/R/examples/austptha_template)
used for the PTHA18.

The study results are provided under a [Creative Commons 4.0 International
Licence](http://creativecommons.org/licenses/by/4.0/legalcode), while the
source-code is provided under a [BSD3 license](../LICENSE). 

Geoscience Australia has tried to make the information in this product as
accurate as possible. However, it does not guarantee that the information is
totally accurate or complete. Therefore, you should not solely rely on this
information when making a commercial decision. 

# Obtaining basic tsunami hazard information.

## Obtaining tsunami maximum-stage exceedance-rates at sites around Australia

The tsunami 'maximum-stage' is the maximum water-level that a particular
tsunami attains at a particular location. This gives an idea of how 'big' the
tsunami is. In the PTHA18 we ignore tidal variations and assume a constant
mean-sea-level (MSL=0), so the 'maximum-stage' is equivalent to the maximum
elevation of the tsunami wave above MSL. 

The maximum-stage exceedance-rates describe 'how often' tsunamis occur with
maximum-stage above a particular threshold value. For example, you could ask how
often maximum-stages above 0.5m (or 2.0m) occur, in terms of the average number of
events each year at a particular site. At most locations there would be less
than one event each year. Thus exceedance-rates are often small numbers.
For example, an exceedance-rate of 0.002=1/500 would correspond to one event
every 500 years on average. Alternatively, one could say the event has a 500
year Average Recurrence Interval (ARI).

In PTHA18 this information is stored at a set of points in the ocean.  These
are termed 'hazard points' because they record site-specific hazard
information. The wave time-series for every scenario can also be obtained at
all hazard points ([see here](DETAILED_README.md)).

Most hazard points are concentrated around Australia and its
territories. We also store DART buoy locations (deep-ocean
gauges which measure tsunami wave heights), which are useful for
model testing. Further, we store a 'thin' layer of hazard points
globally at the 100m contour, which is useful for model testing and inter-comparison.
**If using points far from Austraila, beware that we ignore tsunamigenic
source zones which are not considered relevant for Australia (e.g. in the
Carribbean, the Mediterrean, the Manila trench, Kaikoura in New Zealand,
western Japan). Outside of Australia, you should very carefully consider
how the results can be used, noting they may ignore the most relevant
source-zones.**

The tsunami maximum-stage provides at best a rough characterisation of the
tsunami's capacity to produce coastal inundation. The tsunami inundation will
be affected by the full details of the tsunami wave train, and how it interacts
with the coastal landscape. However, all else being equal, a larger
maximum-stage will generally lead to larger inundation.

The simplest way to examine the PTHA18 tsunami maximum-stage exceedance-rates
is to [download this csv file](http://dapds00.nci.org.au/thredds/fileServer/fj6/PTHA/AustPTHA_1/EVENT_RATES/revised1_tsunami_stages_at_fixed_return_periods.csv).
It contains the following columns:
* `lon`, `lat` give the hazard point location in longitude/latitude (degrees). 
* `elev` is the bathymetry at the hazard point (negative = below MSL)
* `gaugeID` is a hazard point ID (real number).
* multiple columns with names like `STAGE_XXXX` where XXXX is a number, and 1/XXXX is the exceedance-rate. These corresponds to the tsunami maximum-stage which has mean exceedance-rate = 1/XXXX. For example, the column `STAGE_100` gives the tsunami maximum-stage that is exceeded once every 100 years on average, according to the mean of all the rate models in our logic-tree.
* multiple columns with names like `STAGE_upper_ci_XXXX`. These values are similar to the above, but describe the upper limit of the 95% credible interval for the stage with the specified exceedance-rate. (i.e. 97.5 percentile)
* multiple columns with names like `STAGE_lower_ci_XXXX`. These are similar to the above, but describe the lower limit of the 95% credible interval for the stage with the specified exceedance-rate. (i.e. 2.5 percentile)
* multiple columns with names like `STAGE_median_XXXX`. These are similar to the above, but describe the 'epistemic median' stage with the specified exceedance-rate (i.e. 50th percentile)
* multiple columns with names like `STAGE_16pc_XXXX`. These are similar to the above, but describe the 16th percentile.
* multiple columns with names like `STAGE_84pc_XXXX`. These are similar to the above, but describe the 84th percentile.

Note 'max stage' values below 2cm (or above 20m) are treated as missing data
(NA). Such values are unlikely to be of interest, but if necessary they can be
reconstructed from the [detailed information](DETAILED_README.md).

[Similar data is available in shapefile format here](http://dapds00.nci.org.au/thredds/fileServer/fj6/PTHA/AustPTHA_1/EVENT_RATES/revised1_tsunami_stages_at_fixed_return_periods.zip). You must unzip the file after download.
Shapefiles have a weakness; attribute names must not exceed 10 characters. 
Therefore the attributes are renamed in some instances, as compared with the
above csv:
* `lon`, `lat` give the location in longitude/latitude (degrees). 
* `elev` is the bathymetry at the hazard point (negative = below MSL)
* `gaugeID` is a real hazard point ID
* `ST_XXXX` is the same as `STAGE_XXXX` described above
* `STu_XXXX` is the same as `STAGE_upper_ci_XXXX` described above
* `STl_XXXX` is the same as `STAGE_lower_ci_XXXX` described above
* `ST50_XXXX` is the same as `STAGE_median_XXXX` described above
* `ST16_XXXX` is the same as `STAGE_16pc_XXXX` described above
* `ST84_XXXX` is the same as `STAGE_84pc_XXXX` described above

At most hazard points there is large uncertainty in the maximum-stage for a
given exceedance-rate. This is largely due to uncertainty in the frequencies of
high-magnitude subduction zone earthquakes. A more detailed discussion of these
topics can be found in the 
[Australian Tsunami Hazard Modelling Guidelines](https://knowledge.aidr.org.au/media/5640/tsunami-planning-guidelines.pdf).


*Note:* The above results follow the methodology in [this PAGEOPH publication](https://link.springer.com/article/10.1007/s00024-019-02299-w) to compute exceedance-rate percentiles, which is slightly different to the methodology in the original [PTHA18 report](http://dx.doi.org/10.11636/Record.2018.041) (see discussion in Section 3.5 of the PAGEOPH paper). Differences are generally small and unlikely to be important, but for reference the older results can still be obtained in [csv form](http://dapds00.nci.org.au/thredds/fileServer/fj6/PTHA/AustPTHA_1/EVENT_RATES/tsunami_stages_at_fixed_return_periods.csv) and [shapefile form](http://dapds00.nci.org.au/thredds/fileServer/fj6/PTHA/AustPTHA_1/EVENT_RATES/tsunami_stages_at_fixed_return_periods.zip).

## Obtaining site-specific hazard information (including source deaggregation)

For each hazard point, the PTHA18 includes a standard pdf plot which shows:
1. The maximum-stage vs exceedance-rate
2. A convergence check on the above
3. The hazard deaggregation information for a range of return periods
4. Information on the maximum-stage for each unit-source tsunami.

An example plot can be downloaded
[here](http://dapds00.nci.org.au/thredds/fileServer/fj6/PTHA/AustPTHA_1/EVENT_RATES/revised1_station_summary_plots/example_plot.pdf).
Because there are thousands of hazard points, the plots at other sites are provided in a set
of zip folders. Each zip folder contains around 200 sites in a particular
longitude range. The zip folders can be accessed 
[here](http://dap.nci.org.au/thredds/remoteCatalogService?catalog=http://dapds00.nci.org.au/thredds/catalog/fj6/PTHA/AustPTHA_1/EVENT_RATES/revised1_station_summary_plots/catalog.xml).
Follow the link to the http download to get the file. The zip folder names are
of the form *station_summary_plots_longitudes_LOWER_UPPER.zip* where *LOWER* is
the lower longitude limit, and *UPPER* is the upper longitude limit. 

For example if I were searching for a hazard point at (lon,lat)= (151.408,
-34.08), then by inspection of the *LOWER* and *UPPER* longitudes in files
[here](http://dap.nci.org.au/thredds/remoteCatalogService?catalog=http://dapds00.nci.org.au/thredds/catalog/fj6/PTHA/AustPTHA_1/EVENT_RATES/revised1_station_summary_plots/catalog.xml),
it should be contained in the file *station_summary_plots_151.38_152.zip*
(because the *LOWER* and *UPPER* longitudes bracket the value 151.408, which is
the one I want).

*Note:* The exceedance-rate percentile calculation has been revised follow the methodology in [this PAGEOPH publication](https://link.springer.com/article/10.1007/s00024-019-02299-w), as discussed above for the csv and shapefile outputs. The updates are relatively minor, but [the older results are still available here](http://dap.nci.org.au/thredds/remoteCatalogService?catalog=http://dapds00.nci.org.au/thredds/catalog/fj6/PTHA/AustPTHA_1/EVENT_RATES/station_summary_plots/catalog.xml) in case you need to check them for some reason.

## Interpreting exceedance-rate information

Please read the [PTHA18 report](http://dx.doi.org/10.11636/Record.2018.041) and the 
[PAGEOPH Paper](https://link.springer.com/article/10.1007/s00024-019-02299-w) for
further information on interpreting the exceedance-rate information. 

The maximum-stage exceedance-rates vary from site to site, depending on exposure
to earthquake-generated tsunamis. For a given exceedance-rate, there is also a
general tendency for the tsunami size to increase in shallower water. Such
'shoaling' is a well known property of ocean waves. 

The model results are not expected to be accurate everywhere, but **in general
results far offshore and in deep water are expected to be higher quality than
nearshore results**. The reasons are:
* The PTHA18 tsunami model has a spatial grid size of 1 arc minute (around 1.8 km),
and is run on relatively coarse elevation data (a combination of the 
[Australian Bathymetry and Topography Grid 2009](http://www.ga.gov.au/metadata-gateway/metadata/record/gcat_67703)
product, and the global [GEBCO
2014](https://www.gebco.net/data_and_products/gridded_bathymetry_data/) bathymetry grid).
While appropriate for modelling oceanic-scale tsunami propagation, it is not
expected to accurately model tsunamis near the coast and in shallow waters.
* In shallower waters, where wave heights become an appreciable fraction of the water depth, 
the assumptions underlying our linear tsunami model are violated. This is most
likely to be a problem in shallow waters, and for larger tsunamis.

Because of this, **for modelling purposes we strongly encourage the use of
points well offshore in deep water**. If you can, use sites where wave heights of interest
do not exceed a few percent of the water depth. For tsunami propagation modelling,
it may be preferable to simulate the tsunami from source (using initial
conditions [provided here](DETAILED_README.md)), which circumvents these
issues. Modelling from source also also facilitates the use of alternative hydrodynamic
models (e.g. with dispersion and/or friction, which are important in some situations), and 
using other bathymetric data when simulating these large scales. Nearshore
points should **only** be used as a rough guide to possible tsunami wave
heights, **NOT FOR FORCING INUNDATION MODELS**, and should be refined in future
using higher resolution models and data. 

The PTHA18 can help enable national consistency in site-specific tsunami
inundation hazard studies. At any particular coastal site, the tsunami
inundation hazard will reflect both the frequency with which significant
"offshore" tsunamis occur, as well as how these tsunamis interact with the
local coastal topography. The PTHA18 was created to provide a nationally
consistent view of the "offshore" part of the problem. To determine the onshore
hazard at a particular site, the offshore results can be used for force
high-resolution inundation models (which also require detailed site-specific
elevation data). This is analogous to how national scale extreme-rainfall
information provides a consistent basis for flood hazard assessments. 

# Obtaining detailed information on earthquake scenarios, tsunami initial conditions, and wave time-series

Please see the [detailed README](DETAILED_README.md) for information on
extracting this kind of information from the output files.
