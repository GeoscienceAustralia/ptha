# **Guide to accessing the 2018 Australian Probabilistic Tsunami Hazard Assessment (PTHA) results**

# ***NOTE: Currently the study is still in-progress, and results available here will be changed without warning.***

We provide access to basic tsunami hazard information in easy-to-use csv and
shapefile formats. This information is useful to get a high-level overview of
the PTHA results. 

We also provide access to detailed information for every tsunami scenario in our
analysis, including: the earthquake discretization, the tsunami initial
condition, and the resulting wave time-series at a large set of points in the
ocean ('hazard points'). This information is useful for site-specific tsunami
hazard studies (which use it for model 'boundary conditions' and to estimate
how often hazardous tsunamis might occur). Working with this information is
more complex, and is explained in our [detailed README](DETAILED_README.md).

The methodologies used in this study, and associated testing, are discussed in
[the project report](PROVIDE LINK WHEN IT IS PUBLISHED). In addition, the codes
used to conduct the analysis are available open-source in the 
[ptha package](https://github.com/GeoscienceAustralia/ptha).
While that package mostly contains generic functionality for PTHA, there is also
a folder with 
[project specific scripts](https://github.com/GeoscienceAustralia/ptha/tree/master/R/examples/austptha_template)
used for the 2018 Australian PTHA.

The study results are provided under a [Creative Commons 4.0 International
Licence](http://creativecommons.org/licenses/by/4.0/legalcode), while the
source-code is provided under a [BSD3 license](../LICENSE). Geoscience
Australia has tried to make the information in this product as accurate as
possible. However, it does not guarantee that the information is totally
accurate or complete. Therefore, you should not solely rely on this information
when making a commercial decision.


# Obtaining basic tsunami hazard information.

## Obtaining tsunami maximum-stage exceedance-rates at sites around Australia

The tsunami 'maximum-stage' is the maximum water-level that a particular
tsunami attains at a particular location. This gives an idea of how 'big' the
tsunami is. In the current analysis we ignore tidal variations and assume a
constant mean-sea-level (MSL=0), so the 'maximum-stage' is equivalent to the
maximum elevation of the tsunami wave above MSL. 

The maximum-stage exceedance-rates describe 'how often' tsunamis occur with
maximum-stage above a particular threshold value. For example, you could ask how
often maximum-stages above 0.5m (or 2.0m) occur, in terms of the average number of
events each year at a particular site. At most locations there would be less
than one event each year, so the exceedance-rates are typically small numbers.
For example, an exceedance-rate of 0.002=1/500 would correspond to one event every 500
years on average. Alternatively, one could say the event has a 500 year Average
Recurrence Interval (ARI).

In the 2018 PTHA, this information is stored at a set of points in the ocean.
These points are termed 'hazard points' because we provide the hazard
information at these sites. (The wave time-series for every scenario can also be
obtained at all hazard points using methods described later in this document).

Most of the hazard points are concentrated around Australia and its
territories. We also store points at the locations of DART buoys (deep-ocean
gauges which measure tsunami wave heights), as these are very useful for
testing the model. In addition we store a 'thin' layer of hazard points
globally at the 100m contour, which is useful for testing our model and
comparing with previous studies. **If using points far from Austraila, consider
that we ignore tsunamigenic source zones that are not considered relevant for
Australia (e.g. in the Carribbean, the Mediterrean, the Manila trench,
Kaikoura in New Zealand, western Japan). Therefore, outside of Australia, you
should very carefully consider whether the results can be used, noting they may
ignore the most relevant source-zones.**

It should be noted that the tsunami maximum-stage provides only a rough
characterisation of its capacity to cause coastal inundation. Ultimately
tsunami inundation will be affected by the full details of the tsunami wave
train, and how it interacts with the coastal landscape. But all else being
equal, a larger maximum-stage will generally lead to larger inundation.

The simplest way to examine the tsunami maximum-stage exceedance-rates in
the 2018 PTHA is to [download this csv file](http://dapds00.nci.org.au/thredds/fileServer/fj6/PTHA/AustPTHA_1/EVENT_RATES/tsunami_stages_at_fixed_return_periods.csv).
This csv file contains the following columns:

* `lon`, `lat` give the hazard point location in longitude/latitude (degrees). 

* `elev` is the bathymetry at the hazard point (negative = below MSL)

* `gaugeID` is a real hazard point ID

* multiple columns with names like `STAGE_XXXX` where XXXX is a number, and 1/XXXX is the exceedance-rate. These values corresponds to the tsunami maximum-stage which has mean exceedance-rate = 1/XXXX. For example, the column `STAGE_100` gives the tsunami maximum-stage that is exceeded once every 100 years on average, according to the mean of all the rate models in our logic-tree.

* multiple columns with names like `STAGE_upper_ci_XXXX`. These values are similar to the above, but describe the upper limit of the 95% credible interval for the stage with the specified exceedance-rate. (i.e. 97.5 percentile)

* multiple columns with names like `STAGE_lower_ci_XXXX`. These are similar to the above, but describe the lower limit of the 95% credible interval for the stage with the specified exceedance-rate. (i.e. 2.5 percentile)

* multiple columns with names like `STAGE_median_XXXX`. These are similar to the above, but describe the 'epistemic median' stage with the specified exceedance-rate (i.e. 50th percentile)

* multiple columns with names like `STAGE_16pc_XXXX`. These are similar to the above, but describe the 16th percentile.

* multiple columns with names like `STAGE_84pc_XXXX`. These are similar to the above, but describe the 84th percentile.

Note that 'max stage' values below 2cm (or above 20m) are treated as missing
data in this file (NA). While such values are unlikely to be of interest, if
necessary they can be reconstructed from the detailed information we provide
(later in this document).

[Similar data is available in shapefile format here](http://dapds00.nci.org.au/thredds/fileServer/fj6/PTHA/AustPTHA_1/EVENT_RATES/tsunami_stages_at_fixed_return_periods.zip). You will need to unzip the file after download.
A shortcoming of the shapefile format is that there is a 10 character limit on
attribute names. Therefore the attributes are renamed in some instances, as
compared with the above csv:

* `lon`, `lat` give the location in longitude/latitude (degrees). 

* `elev` is the bathymetry at the hazard point (negative = below MSL)

* `gaugeID` is a real hazard point ID

* `ST_XXXX` is the same as `STAGE_XXXX` described above

* `STu_XXXX` is the same as `STAGE_upper_ci_XXXX` described above

* `STl_XXXX` is the same as `STAGE_lower_ci_XXXX` described above

* `ST50_XXXX` is the same as `STAGE_median_XXXX` described above

* `ST16_XXXX` is the same as `STAGE_16pc_XXXX` described above

* `ST84_XXXX` is the same as `STAGE_84pc_XXXX` described above

At most hazard points you will find there is large uncertainty in the
maximum-stage for a given exceedance-rate. This is mainly due to large uncertainty
in the frequencies of high-magnitude subduction zone earthquakes. A more
detailed discussion of these topics can be found in the 
[Australian Tsunami Hazard Modelling Guidelines](https://knowledge.aidr.org.au/media/5640/tsunami-planning-guidelines.pdf).

## Obtaining site-specific hazard information (including source deaggregation)

For each hazard point we provide a standard pdf plot which shows:
1. The maximum-stage vs exceedance-rate
2. A convergence check on the above
3. The hazard deaggregation information for a range of return periods
4. Information on the maximum-stage for each unit-source tsunami.

An example plot can be downloaded
[here](http://dapds00.nci.org.au/thredds/fileServer/fj6/PTHA/AustPTHA_1/EVENT_RATES/station_summary_plots/example_plot.pdf).
Because there are thousands of hazard points, the plots at other sites are provided in a set
of zip folders. Each zip folder containing around 200 sites in a particular longitude
range. The zip folders can be accessed 
[here](http://dap.nci.org.au/thredds/remoteCatalogService?catalog=http://dapds00.nci.org.au/thredds/catalog/fj6/PTHA/AustPTHA_1/EVENT_RATES/station_summary_plots/catalog.xml).
Follow the link to the http download to get the file. The zip folder names are
of the form *station_summary_plots_longitudes_LOWER_UPPER.zip* where *LOWER* is
the lower longitude limit, and *UPPER* is the upper longitude limit. 

For example if I were searching for a hazard point at (lon,lat)= (151.408,
-34.08), then by inspection of the *LOWER* and *UPPER* longitudes in files
[here](http://dap.nci.org.au/thredds/remoteCatalogService?catalog=http://dapds00.nci.org.au/thredds/catalog/fj6/PTHA/AustPTHA_1/EVENT_RATES/station_summary_plots/catalog.xml),
it should be contained in the file *station_summary_plots_151.38_152.zip*
(because the *LOWER* and *UPPER* longitudes bracket the value 151.408, which is
the one I want).

## Interpreting exceedance-rate information

Please read the [PTHA18 report](PROVIDE LINK WHEN AVAILABLE) for further
information on interpreting the exceedance rate information. 

The maximum-stage exceedance-rates vary from site to site, depending on exposure
to earthquake-generated tsunamis. For a given exceedance-rate, there is also a
general tendency for the tsunami size to increase in shallower water. Such
'shoaling' is a well known property of ocean waves. 

The model results are not expected to be accurate everywhere, but **in general
results far offshore and in deep water are expected to be higher quality than
nearshore results**. The reasons are:

* Our tsunami model has a spatial grid size of 1 arc minute (around 1.8 km), 
and is run on relatively coarse elevation data (a combination of the 
[Australian Bathymetry and Topography Grid 2009](http://www.ga.gov.au/metadata-gateway/metadata/record/gcat_67703)
product, and the global [GEBCO
2014](https://www.gebco.net/data_and_products/gridded_bathymetry_data/) bathymetry grid).
While appropriate for modelling oceanic-scale tsunami propagation, it is not
expected to accurately model tsunamis near the coast and in shallow waters.

* At locations where wave heights become an appreciable fraction of the water depth, 
the modelled waves will violate the assumptions underlying our linear tsunami
model. This is most likely to be a problem in shallow waters.

Because of this, **for modelling purposes we strongly encourage the use of
points well offshore in deep water** (preferably with wave heights of interest
not exceeding a few percent of the water depth). For tsunami propagation modelling,
we encourage users to simulate the from source using the initial conditions we provide,
which circumvents these issues. Nearshore points should only be used as a rough
guide to possible tsunami wave heights, and should be refined in future using
higher resolution models and data. 

The above discussion might lead non-specialists to wonder why we develop the
PTHA at offshore points - considering tsunami inundation is of most interest
for risk management. The key motivation is to provide a nationally consistent
source of tsunami scenarios and exceedance-rates for local-scale tsunami hazard
studies. At any particular coastal site, the tsunami inundation hazard will
reflect both the frequency with which significant "offshore" tsunamis occur, as
well as how these tsunamis interact with the local coastal topography. The PTHA
enables a nationally consistent view of the "offshore" parts of the problem. On
the other hand, site-by-site assessments are generally required to determine
the interactions of tsunamis with the coastal topography (mainly because it
involves computationally intensive modelling with detailed input data). 

In this way, the PTHA facilitates national consistency in site-specific tsunami
inundation hazard studies. This is analogous to how national scale
extreme-rainfall information provides a consistent basis for flood hazard
assessments. 

# Obtaining detailed information on earthquake scenarios, tsunami initial conditions, and wave time-series

Please see the [detailed README](DETAILED_README.md) for information on
extracting this kind of information from the output files.
