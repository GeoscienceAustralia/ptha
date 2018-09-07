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
how often hazardous tsunamis might occur). To access this information users
will need to install a range of software which allows working with the PTHA
data on the NCI THREDDS server. 

The methodologies used in this study, and associated testing, are discussed in
the project report [provide link when it is published][]. In addition, the codes
used to conduct the analysis are available open-source in the 
[ptha package](https://github.com/GeoscienceAustralia/ptha).
While that package mostly contains generic functionality for PTHA, there is also
a folder with 
[project specific scripts](https://github.com/GeoscienceAustralia/ptha/tree/master/R/examples/austptha_template)
used for the 2018 Australian PTHA.

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

* multiple columns with names like `STAGE_upper_ci_XXXX`. These values are similar to the above, but describe the upper limit of the 95% credible interval for the stage with the specified exceedance-rate. (i.e. 97.5% quantile)

* multiple columns with names like `STAGE_lower_ci_XXXX`. These are similar to the above, but describe the lower limit of the 95% credible interval for the stage with the specified exceedance-rate. (i.e. 2.5% quantile)

* multiple columns with names like `STAGE_median_XXXX`. These are similar to the above, but describe the 'epistemic median' stage with the specified exceedance-rate (i.e. 50% quantile)

* multiple columns with names like `STAGE_16pc_XXXX`. These are similar to the above, but describe the 16% quantile.

* multiple columns with names like `STAGE_84pc_XXXX`. These are similar to the above, but describe the 84% quantile.

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

## Obtaining more detailed exceedance-rate information for specific sites

FIXME describe where to obtain the pdf files containing the summary information

## Interpreting exceedance-rate information

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
not exceeding a few percent of the water depth). Nearshore points should only
be used as a rough guide to possible tsunami wave heights, and should be
refined in future using higher resolution models and data. 

The above statements might lead non-specialists to question the purpose of this
PTHA, given that for risk management purposes the tsunami inundation is of most
interest. The key reason for developing an 'offshore' PTHA is that it provides
essential input data to support the high-resolution inundation models required
for tsunami risk management. The tsunami scenarios and associated
exceedance-rates provided in this PTHA will be used as 'boundary conditions' to
drive the site-specific tsunami inundation models. This helps facilitate
national consistency in tsunami inundation modelling, while reducing the need
for nearshore tsunami modellers to develop expertise in subjects such as
earthquake kinematics and tsunami generation, earthquake magnitude-frequency
relations, and quantification of the associated uncertainties. Furthermore, we
have been able to extensively test the performance of our tsunami generation
methods using global-scale deep ocean tsunami observations for 20 historical
tsunami events. Most local-scale studies would not have the time or resources
to undertake such testing as part of their project.

# Obtaining detailed information on earthquake scenarios, tsunami initial conditions, and wave time-series

For every scenario in our analysis we provide earthquake information, tsunami
initial conditions, and wave time-series at every hazard point. Combined with
the exceedance-rate modelling, such inputs can be used to drive local scale
tsunami inundation models for hazard and risk assessments.

To access the detailed information, the user needs to interact with our files
via the NCI THREDDS server. We provide R scripts to facilitate this, and the
process is described below. A range of software must be installed to run these
codes, [as described here](INSTALL.md)

Unfortunately the installation and data extraction process may be challenging
for users with limited experience in scientific programming and Linux. Users
doing tsunami hazard studies **in Australia** can alternatively contact
Geoscience Australia directly if they have difficulty with any of these steps
(please email Gareth Davies at gareth.davies@ga.gov.au). 

## **Usage**

Make sure you have successfully installed the software [as described here](INSTALL.md).
Please confirm that everything installed correctly by running the [test_all.R](test_all.R) script.
You need to start R in the directory that contains that script and this README
file (i.e. the directory name will end in `/ptha/ptha_access/`)

```r
# This should print 'PASS' a few times. If not, something is wrong with your
# install (or perhaps your internet connection is slow -- give it a few tries).
source('test_all.R')
```

```
## [1] "PASS"
## [1] "PASS"
## [1] "PASS"
## [1] "PASS"
## [1] "PASS"
## [1] "PASS"
## [1] "PASS"
## [1] "PASS"
```

### ***Viewing the locations of hazard points and source zones***

It is possible to view the hazard points from an interactive map in R. Note however
that similar information is provided in the csv and shapefiles above, and for
many users it will be easier to view those files using GIS and/or a spreadsheet
application.

To view the source-zones and hazard points on an interactive map in R, start
R in the same directory that the [hazard_points_plot.R](hazard_points_plot.R)
file resides in, and do:

```r
source('hazard_points_plot.R')
```
The should open a map in your web browser, containing all unit sources and
hazard points. 

The first time you run this code it will download several datasets to your
machine for use in the map. These will be placed in the DATA and SOURCE_ZONES
folders. The download might take a minute or more, depending on your internet
connection. Future runs will read the data from your machine, so should be
faster. If you want to download fresh data (e.g. if you think there has been an
update, or your files seem corrupted), then just manually delete the DATA
and SOURCE_ZONES folders.


![hazardpoints1](figure/hazard_point_viewer_screenshot1.png)

Initially, most of the hazard points will be aggregated into coloured circles
containing clusters of hazard points. This is done because it is too slow to
render all hazard points at the same time on the one map. In the above figure,
we see green circles (containing less than 10 hazard points), yellow circles
(containing 10-100 hazard points), and red circles (containing more than 100
hazard points). A number on the circle shows how many hazard points they
contain. There are also a few individual hazard points (which are far from
others), and in the above figure they mostly correspond to the locations of DART
buoys.

If you zoom in enough (e.g. below we look at Christmas Island), eventually the circles
containing many points should be replaced by individual hazard points
(circles). They can be queried with a mouse click. For each point, we store
basic stage-vs-exceedance-rate information, as was discussed above. Note stage values
below 2cm or above 20m are reported as NA.
![hazardpoints2](figure/hazard_point_viewer_screenshot2.png)

The unit sources appear as a polygonal grid. Individual unit sources can also
be queried (e.g. to learn the name of the source-zone in our analysis) 
![hazardpoints3](figure/hazard_point_viewer_screenshot3c.png)
The controls on the top left of the map can be expanded as shown in the figure.
These should allow you to change the background layer, and to turn layers on
and off.

### ***Obtaining metadata on the earthquake scenarios on each source-zone***

Earthquake scenario metadata is accessed on a per-source-zone basis. In a typical
application you would use the detailed exceedance-rate plots discussed above
to identify the main source-zones of interest for a particular site. Below
we show an example using the `puysegur` source-zone, which is located just south
of New Zealand, to the north of Macquarie Island.

To download metadata from the NCI describing the earthquake scenarios on a
particular source-zone, start R in the current directory, and do:

```r
# Import the functions
source('get_PTHA_results.R')

# Find the possible names of the source-zones
get_source_zone_events_data()
```

```
## [1] "You did not pass a valid source_zone to get_source_zone_events_data. The allowed source_zone values are:"
##  [1] "   alaskaaleutians"         "   arutrough"              
##  [3] "   banda_detachment"        "   cascadia"               
##  [5] "   floreswetar"             "   hjort"                  
##  [7] "   izumariana"              "   kermadectonga2"         
##  [9] "   kurilsjapan"             "   macquarieislandnorth"   
## [11] "   makran2"                 "   manokwari"              
## [13] "   manus"                   "   mexico"                 
## [15] "   moresby_trough"          "   mussau"                 
## [17] "   newguinea2"              "   newhebrides2"           
## [19] "   north_sulawesi"          "   outer_rise_timor"       
## [21] "   outerrise_kermadectonga" "   outerrise_puysegur"     
## [23] "   outerrisenewhebrides"    "   outerrisesolomon"       
## [25] "   outerrisesunda"          "   philippine"             
## [27] "   puysegur2"               "   ryuku"                  
## [29] "   sandwich"                "   sangihe"                
## [31] "   sangihe_backthrust"      "   se_sulawesi"            
## [33] "   seram_thrust"            "   seramsouth"             
## [35] "   solomon2"                "   southamerica"           
## [37] "   sunda2"                  "   tanimbar"               
## [39] "   timortrough"             "   tolo_thrust"            
## [41] "   trobriand"              
## [1] "Please pass one of the above source_zone names to this function to get its metadata"
```

Above we called the function `get_source_zone_events_data` that is typically
used to get metadata on all scenarios on a particular source-zone. However, if
no arguments are passed, then by default it prints the valid `source_zone` names
and exits. That can help you learn what the source-zone names are. They can
also be found by clicking on source-zones in the interactive hazard map
(discussed above).

Suppose we are interested in the Puysegur source-zone. From the above list and/or the 
interactive viewer, we could infer that it was named `puysegur2`. We can then get the
scenario metadata as follows:

```r
# Example: get metadata for the puysegur source_zone
puysegur = get_source_zone_events_data('puysegur2')
```

This variable `puysegur` is now an R `list`, which contains two `data.frame`'s
summarising the source-zone geometry and the earthquake secnarios, and a character
vector giving the associated tide-gauge files (where tsunami time-series are
stored)

```r
names(puysegur)
```

```
## [1] "events"                 "unit_source_statistics"
## [3] "gauge_netcdf_files"
```

```r
lapply(puysegur, class) # Get class of each entry in the list 'puysegur'
```

```
## $events
## [1] "data.frame"
## 
## $unit_source_statistics
## [1] "data.frame"
## 
## $gauge_netcdf_files
## [1] "character"
```

We now describe the unit-source-statistics table.
`puysegur$unit_source_statistics` contains summary statistics about the
unit-sources. For each unit source this gives the centroid `lon` and `lat` and
`depth`; the unit source dimensions `length` and `width`; the rupture source
mechanism (`strike`, `dip`, `rake`); and indices `downdip_number`,
`alongstrike_number`, and `subfault_number` which give information of the
placement of the unit source on the grid of all unit sources.

```r
# Get the names of all summary statistics
names(puysegur$unit_source_statistics)
```

```
##  [1] "lon_c"                  "lat_c"                 
##  [3] "depth"                  "strike"                
##  [5] "dip"                    "rake"                  
##  [7] "slip"                   "length"                
##  [9] "width"                  "downdip_number"        
## [11] "alongstrike_number"     "subfault_number"       
## [13] "max_depth"              "initial_condition_file"
## [15] "tide_gauge_file"
```
Here we determine the dimensions of the table, and look at a few rows

```r
# Get the table dimensions
dim(puysegur$unit_source_statistics)
```

```
## [1] 28 15
```

```r
# Print rows 1 and 2
puysegur$unit_source_statistics[1:2,]
```

```
##      lon_c     lat_c     depth   strike      dip rake slip   length
## 1 163.6647 -49.88909  6.364741 20.39946 18.67347   90    1 45.92458
## 2 164.1119 -49.96848 26.364741 14.18287 43.24057   90    1 45.54654
##      width downdip_number alongstrike_number subfault_number max_depth
## 1 39.80127              1                  1               1  12.83707
## 2 39.98489              2                  1               2  40.00000
##                                                                                        initial_condition_file
## 1 /g/data1a/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/puysegur2/EQ_SOURCE/Unit_source_data/puysegur2/puysegur2_1_1.tif
## 2 /g/data1a/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/puysegur2/EQ_SOURCE/Unit_source_data/puysegur2/puysegur2_2_1.tif
##                                                                                                                                                                         tide_gauge_file
## 1 /g/data/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/puysegur2/TSUNAMI_UNIT_SOURCES/unit_source_tsunami/RUN_20180815145902_puysegur2_1_1/RUN_ID100001_20180816_044751.819/Gauges_data_ID100001.nc
## 2 /g/data/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/puysegur2/TSUNAMI_UNIT_SOURCES/unit_source_tsunami/RUN_20180815145906_puysegur2_2_1/RUN_ID100001_20180816_050648.956/Gauges_data_ID100001.nc
```
In addition, the `initial_condition_file` and `tide_gauge_file` variables
provide a link to the vertical deformation and tsunami model run respectively,
for each unit source. Note that the file path names sometimes differ slightly
from the location on the NCI Thredds server (although in a fairly obvious way,
e.g. sometimes we see /g/data replaced with /g/data1a). The scripts provided
here will translate the file paths as required for remote access.

Next we consider the scenario metadata table. `puysegur$events` contains summary
statistics about the earthquake scenarios. 

```r
# Print the names of all scenario summary statistics
names(puysegur$events)
```

```
##  [1] "event_index_string"                  
##  [2] "event_slip_string"                   
##  [3] "Mw"                                  
##  [4] "target_lon"                          
##  [5] "target_lat"                          
##  [6] "peak_slip_downdip_ind"               
##  [7] "peak_slip_alongstrike_ind"           
##  [8] "physical_corner_wavenumber_x"        
##  [9] "physical_corner_wavenumber_y"        
## [10] "sourcename"                          
## [11] "uniform_event_row"                   
## [12] "rate_annual"                         
## [13] "rate_annual_lower_ci"                
## [14] "rate_annual_upper_ci"                
## [15] "variable_mu_Mw"                      
## [16] "variable_mu_rate_annual"             
## [17] "variable_mu_rate_annual_lower_ci"    
## [18] "variable_mu_rate_annual_upper_ci"    
## [19] "variable_mu_rate_annual_median"      
## [20] "variable_mu_rate_annual_16pc"        
## [21] "variable_mu_rate_annual_84pc"        
## [22] "variable_mu_weight_with_nonzero_rate"
## [23] "weight_with_nonzero_rate"            
## [24] "rate_annual_16pc"                    
## [25] "rate_annual_84pc"                    
## [26] "rate_annual_median"
```

```r
# Get the table dimensions
dim(puysegur$events)
```

```
## [1] 5881   26
```
While there are many ways to investigate the `events` table, a simple approach is
to just print some rows. In general low row-indices will correspond to low
magnitudes, and high indices to high magnitudes.


```r
# Print some rows (we choose 2050, 2051, 2052)
puysegur$events[2050:2052, ]
```

```
##      event_index_string        event_slip_string  Mw target_lon target_lat
## 2050       21-23-25-27- 5.838_8.941_2.336_1.508_ 7.9   166.7325  -45.56029
## 2051          22-23-24-        4.554_5.927_8.27_ 7.9   166.7325  -45.56029
## 2052       19-20-21-22-   1.42_3.54_2.334_7.838_ 7.9   166.7325  -45.56029
##      peak_slip_downdip_ind peak_slip_alongstrike_ind
## 2050                     1                        12
## 2051                     2                        12
## 2052                     2                        11
##      physical_corner_wavenumber_x physical_corner_wavenumber_y sourcename
## 2050                  0.003639969                  0.004682143  puysegur2
## 2051                  0.009215376                  0.011996502  puysegur2
## 2052                  0.006078464                  0.010341691  puysegur2
##      uniform_event_row  rate_annual rate_annual_lower_ci
## 2050               131 2.699053e-05         1.247440e-06
## 2051               131 2.825624e-05         1.305938e-06
## 2052               131 3.113489e-05         1.438982e-06
##      rate_annual_upper_ci variable_mu_Mw variable_mu_rate_annual
## 2050         4.415802e-05       7.612915            3.805892e-05
## 2051         4.622879e-05       7.977218            3.680982e-05
## 2052         5.093842e-05       7.980054            3.509007e-05
##      variable_mu_rate_annual_lower_ci variable_mu_rate_annual_upper_ci
## 2050                     1.977616e-06                     6.544703e-05
## 2051                     1.912711e-06                     6.329906e-05
## 2052                     1.823349e-06                     6.034173e-05
##      variable_mu_rate_annual_median variable_mu_rate_annual_16pc
## 2050                   3.310225e-05                 4.604935e-05
## 2051                   3.201583e-05                 4.453801e-05
## 2052                   3.052005e-05                 4.245719e-05
##      variable_mu_rate_annual_84pc variable_mu_weight_with_nonzero_rate
## 2050                 4.295668e-05                            0.9956551
## 2051                 4.154684e-05                            0.9956551
## 2052                 3.960577e-05                            0.9956551
##      weight_with_nonzero_rate rate_annual_16pc rate_annual_84pc
## 2050                0.9952935     3.260319e-05     3.074646e-05
## 2051                0.9952935     3.413211e-05     3.218830e-05
## 2052                0.9952935     3.760937e-05     3.546753e-05
##      rate_annual_median
## 2050       2.427181e-05
## 2051       2.541002e-05
## 2052       2.799871e-05
```

The most important variables from a users perspective are the moment magnitude
`Mw`, and the "variable shear modulus" moment magnitude `variable_mu_Mw`.
You may be surprised to see we store two different earthquake magnitudes for
each scenario. The `Mw` column holds the earthquake moment magnitude, derived
under the assumption that the shear modulus (or rigidity, i.e. a material
property of the earth) on the source-zone is constant, with a value of 30 GPa
on thrust sources and 60 GPa on normal sources. These values are quite typical
and produce tsunamis that compare well with our DART buoy test set (20
historical tsunamis). However, on subduction zones there is some evidence that
the shear modulus increases with depth, with particularly low values possible
at shallow depths (which may partially account for so-called 'tsunami
earthquakes', which generate large tsunamis compared to their earthquake
magnitude). To account for this we use a depth varying shear modulus model, and
re-compute the magnitude for each scenario (stored in `variable_mu_Mw`), without
changing any other properties of the earthquake. The effect is that shallow earthquakes
have `variable_mu_Mw` \< `Mw`, while the opposite occurs for deep earthquakes.
Our tsunami scenarios also compare well with the DART buoy dataset using this
redefined magnitude (see the PTHA report). 

Some other important variables are the `event_slip_string` and the
`event_index_string`. These variables can be used to determine which
unit-sources are included in the earthquake, and how much slip they have. Note
they are stored as strings with a separator, to permit efficient storage of
earthquakes with a range of sizes. The integer values in `event_index_string` 
correspond to the `subfault_number` values in the unit-source statistics table
discussed above.

Another useful variable is the `weight_with_nonzero_rate`. This gives the
fraction of the exceedance-rate models in the logic tree that suggest scenarios
with the given `Mw` could possibly occur. Values close to 1.0 indicate "a high
fraction of our rate models suggest scenarios with this `Mw` could occur, given a
long enough time-frame". On the other hand, values close to 0.0 indicate that
"a high fraction of our rate models suggest scenarios with this `Mw` would never
occur", with zero corresponding to an impossible scenario (i.e.  according to the
model).

All of our source-zone `events` tables contain scenarios with `Mw` values
ranging from 7.2 to 9.8. This is done for computational convenience,
irrespective of whether we consider the high magnitude scenarios are possible
on the source-zone.  You will notice that scenarios at very large `Mw` always
have a `weight_with_nonzero_rate` equal to zero.

#### ***Understanding the scenario rates***

There are also a number of variables with names including `_rate_annual`.
**Beware: These do not give the exceedance rates for the scenarios!**. Instead
they represent a (fairly nominal) scenario-specific portion of the source-zone's
magnitude-vs-exceedance-rate curve, evaluated at the logic-tree mean and
various quantiles. 

Assigning rates to scenarios like this turns out to be very useful to
facilitate other calculations. For instance, by adjusting these rates we can
make scenarios more or less likely on some parts of the source-zone (e.g to
reflect spatial variations in tectonic convergence rates).  Furthermore, we can
re-weight the scenarios based on their slip, which is used in the PTHA18 to
adjust for model biases. This is done differently for models for constant and
variable shear modulus. Notice that some of the rate variables begin with
`variable_mu_` (they are the variable shear modulus versions, where the shear
modulus varies with depth) -- while the other ones assume constant shear
modulus. 

To understand the sceanrio rates, the key idea is that **if you sum all of the
scenario rates above a given magnitude, then the result will correspond to the
magnitude-vs-exceedance-rate curve for the source-zone**. For example, to get
the rate of scenarios above magnitude 7.85 on this source-zone, with various
logic-tree percentiles describing the uncertainty, you could do the following
calculations:

```r
# Rate of events with Mw > 7.85 -- logic-tree mean
sum(puysegur$events$rate_annual * (puysegur$events$Mw > 7.85))
```

```
## [1] 0.01089304
```
You can do a similar calculation to recover the percentile curves. We store 
the 2.5% percentile (`_lower`)

```r
# Rate of events with Mw > 7.85 -- logic-tree 2.5 percentile
sum(puysegur$events$rate_annual_lower_ci * (puysegur$events$Mw > 7.85))
```

```
## [1] 0.0001980731
```
... and the 16th percentile (`_16pc`)

```r
# Rate of events with Mw > 7.85 -- logic-tree 16 percentile
sum(puysegur$events$rate_annual_16pc * (puysegur$events$Mw > 7.85))
```

```
## [1] 0.005176857
```
... and the median (`_median`)

```r
# Rate of events with Mw > 7.85 -- logic-tree median
sum(puysegur$events$rate_annual_median * (puysegur$events$Mw > 7.85))
```

```
## [1] 0.01059198
```
... and the 84th percentile (`_84pc`)

```r
# Rate of events with Mw > 7.85 -- logic-tree 84 percentile
sum(puysegur$events$rate_annual_84pc * (puysegur$events$Mw > 7.85))
```

```
## [1] 0.01599209
```
... and the 97.5 percentile (`_upper`)

```r
# Rate of events with Mw > 7.85 -- logic-tree 97.5 percentile
sum(puysegur$events$rate_annual_upper * (puysegur$events$Mw > 7.85))
```

```
## [1] 0.02427671
```
Notice that these are ordered as expected, i.e. 2.5% <= 16% <= median <= 84% <=
97.5%. However, please note that at the individual scenario level this
ordering may not hold. For example, the `rate_annual_16pc` value is sometimes
larger than the `rate_annual_median` value. This is because the latter
variables are related to the *derivatives* of the Mw-exceedance-rate
percentile curves. The derivatives might not be ordered in the same way as the 
exceedance-rate curves themselves.

If you do similar calculations using the `variable_mu` versions of these
variables, you will find the results are not identical, although they are
usually quite similar. This is because shear modulus variability changes the
relationship between earthquake magnitude and average slip, which in turn
affects the relationship between the rate of earthquakes and the implied
tectonic-plate motion rates. See the PTHA18 report for a full explanation of
these issues.

In the previous calculations we used 7.85 as the magnitude threshold. It is a
good idea to avoid using values that correspond to the scenario magnitude, i.e.
`7.2, 7.3, ... 9.8`.  The reason is that the earthquake magnitudes are stored
as floating point values, so there will be some minor round-off in the
magnitudes we read from file. This could influence the results of `>` or `>=`
type operations. To avoid that, it is better to use intermdiate magnitudes like
`7.15, 7.25, ... 9.85`.

### ***Obtaining tsunami initial conditions for a single earthquake-tsunami scenario***

Suppose we want to get the tsunami initial conditions (i.e. water surface
deformation) for the earthquake event on row 2050 of `puysegur$events`. The
metadata for event 2050 is:

```r
row_index = 2050 # Use this variable to refer to event 2050
puysegur$events[row_index,]
```

```
##      event_index_string        event_slip_string  Mw target_lon target_lat
## 2050       21-23-25-27- 5.838_8.941_2.336_1.508_ 7.9   166.7325  -45.56029
##      peak_slip_downdip_ind peak_slip_alongstrike_ind
## 2050                     1                        12
##      physical_corner_wavenumber_x physical_corner_wavenumber_y sourcename
## 2050                  0.003639969                  0.004682143  puysegur2
##      uniform_event_row  rate_annual rate_annual_lower_ci
## 2050               131 2.699053e-05          1.24744e-06
##      rate_annual_upper_ci variable_mu_Mw variable_mu_rate_annual
## 2050         4.415802e-05       7.612915            3.805892e-05
##      variable_mu_rate_annual_lower_ci variable_mu_rate_annual_upper_ci
## 2050                     1.977616e-06                     6.544703e-05
##      variable_mu_rate_annual_median variable_mu_rate_annual_16pc
## 2050                   3.310225e-05                 4.604935e-05
##      variable_mu_rate_annual_84pc variable_mu_weight_with_nonzero_rate
## 2050                 4.295668e-05                            0.9956551
##      weight_with_nonzero_rate rate_annual_16pc rate_annual_84pc
## 2050                0.9952935     3.260319e-05     3.074646e-05
##      rate_annual_median
## 2050       2.427181e-05
```
To get its initial condition, you pass the earthquake metadata to the function
`get_initial_condition_for_event`:

```r
# Get the initial condition as a geo-referenced raster
initial_condition = get_initial_condition_for_event(puysegur, row_index)

## The raster can be save as a geotif for use in other software, with:
# writeRaster(initial_conditions, 'my_output_filename.tif')

# Make a plot
plot(initial_condition, main='Initial water surface deformation \n for the example event, Puysegur')
```

![plot of chunk raster_eventXXX](figure/raster_eventXXX-1.png)

The function `get_initial_condition_for_event` used above will download the
unit-source tsunami deformations included in the event and save them in the
folder `SOURCE_ZONES/puysegur/EQ_SOURCE/Unit_source_data/puysegur`. Then it
will sum them, scaled by the unit-source slip, to produce the initial
condition. Next time you call the function, it will check whether the required
files exist in the local folder, and only download those that it needs.
However, you can force the function to download the files (and overwrite any
existing ones) by adding the argument `force_file_download=TRUE` (by default
the latter is `FALSE`). This is useful if the NCI analysis has been updated, or
if you suspect your files have been corrupted somehow (although we have not
seen that).

```r
# Get the initial condition as a geo-referenced raster, forcing download of
# all files from NCI irrespective of whether they exist on the current
# machine
initial_condition = get_initial_condition_for_event(puysegur, row_index, force_file_download=TRUE)
```


### ***Obtaining hazard curves at a particular hazard point***

FIXME: Integrate with above discussion. Consider showing how to download numeric
curve values for a particular point.


### ***Finding earthquake events within a particular wave-height range at a particular hazard point***

FIXME: 

### ***Extracting the tsunami time-series for a particular event at a particular hazard point***

Here we show how to read a flow time-series for a given earthquake event, at a
given hazard point. To do this, you have to know the hazard point `gaugeID`,
which can be found by examining the maximum-stage vs exceedance-rate datasets (csv
and shapefile), or by using the interactive hazard point viewer above. (***In the
latter case, please do not confuse this with the Feature ID that is shown by
default in the interactive map - I would like to remove this field, but do not
yet know how/if it can be done!***). 


```r
# Get stage, uh, vh time-series at DART gauges 55015 and 55042
# To find the ID's, look on the interactive hazard-point map.
model_ts = get_flow_time_series_at_hazard_point(puysegur, 
    event_ID=row_index, 
    hazard_point_ID=c(55015.4, 55042.4))
# Should have a 'time' vector, and 'flow' list, and a 'locations' data.frame, as
# well as the 'events' data
names(model_ts)
```

```
## [1] "time"      "flow"      "locations" "events"
```

```r
# The 'flow' list should have one matrix for each gauge. 
names(model_ts$flow)
```

```
## [1] "55015.4" "55042.4"
```

```r
# Alternatively the user can keep 'flow' as an array with the first dimension
# size equal to the number of gauges, by passing the argument 'unpack_to_list=FALSE'
# The latter option may be more efficient for some computations.

# By default for each gauge, model_ts$flow[["gauge_id"]] is a 3D array. 
# The first dimension is always length 1, the second dimension has length
# equal to the number of time-steps, and the third dimension is of length
# three -- with 1 = Stage, 2 = UH, 3 = VH
dim(model_ts$flow[['55015.4']])
```

```
## [1]    1 4321    3
```

```r
# Example plot of stage
plot(model_ts$time, model_ts$flow[['55015.4']][1,,1], t='l', 
    xlim=c(0,10000), xlab='Seconds after earthquake', ylab='Stage (m)',
    ylim=c(-0.2, 0.35))
points(model_ts$time, model_ts$flow[['55042.4']][1,,1], t='l', 
    col='red')
legend('topright', c('55015.4', '55042.4'), col=c('black', 'red'), lty=c(1,1))

title('Stage time-series for the event at 2 gauges')
```

![plot of chunk getflow](figure/getflow-1.png)

To export the flow time-series to a csv, you can do something like this for
the station of interest:

```r
# Name the site
sitename = '55015.4'
# Note you can get a vector with all names using the comment:
#    names(model_ts$flow)
# and this will allow programatically working with the names

# Make a data.frame with the required data
site_flow = data.frame(
    time=model_ts$time, 
    stage = model_ts$flow[[sitename]][1,,1],
    uh = model_ts$flow[[sitename]][1,,2],
    vh = model_ts$flow[[sitename]][1,,3])

# Save it to a csv
output_file = paste0('output_gauge_data_puysegur_event_', row_index, '_station_', 
    sitename, '.csv')
write.csv(site_flow, output_file, row.names=FALSE)
```
