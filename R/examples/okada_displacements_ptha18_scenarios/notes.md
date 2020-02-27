
# 1) How do I read the data created by running "get_displacements_for_events.R"?

From within R, do this.


```r
    library(rptha) # Get the package
    load('3D_displacements_R_image.Rdata') # Get the data
```

The above command will make your R session contain a bunch of variables, just
as though you had ran the script yourself. To print the names of these variables, do:

```r
    ls() # Like dir() in python
```

```
##  [1] "clean_gaugeID"                                  
##  [2] "config_env"                                     
##  [3] "easting_disp_files"                             
##  [4] "event_inds_list"                                
##  [5] "event_slip"                                     
##  [6] "event_slips_list"                               
##  [7] "event_unit_source_inds"                         
##  [8] "get_displacements_due_to_unit_source"           
##  [9] "get_flow_time_series_at_hazard_point"           
## [10] "get_flow_time_series_SWALS"                     
## [11] "get_initial_condition_for_event"                
## [12] "get_netcdf_attribute_initial_stage_raster"      
## [13] "get_netcdf_gauge_index_matching_ID"             
## [14] "get_netcdf_gauge_indices_in_polygon"            
## [15] "get_netcdf_gauge_indices_near_points"           
## [16] "get_netcdf_gauge_locations"                     
## [17] "get_netcdf_gauge_output_times"                  
## [18] "get_peak_stage_at_point_for_each_event"         
## [19] "get_source_zone_events_data"                    
## [20] "get_stage_exceedance_rate_curve_at_hazard_point"
## [21] "get_stage_exceedance_rate_curves_all_sources"   
## [22] "get_supporting_data"                            
## [23] "i"                                              
## [24] "kt2"                                            
## [25] "make_tsunami_event_from_unit_sources"           
## [26] "northing_disp_files"                            
## [27] "parse_ID_point_index_to_index"                  
## [28] "possible_inds"                                  
## [29] "sort_tide_gauge_files_by_unit_source_table"     
## [30] "source_zone"                                    
## [31] "summarise_events"                               
## [32] "target_pt"                                      
## [33] "test_sum_tsunami_unit_sources"                  
## [34] "vertical_disp_files"                            
## [35] "xyz_displacement_events"                        
## [36] "xyz_displacement_unit_sources"
```

# 2) How do I associate each 3D displacement vector with the correct earthquake slip distribution?

After reading the following comments, it might also help to read the script [get_displacements_for_events.R](get_displacements_for_events.R).
  
## Getting the earthquake events and their slip distributions

The earthquake events are stored in a list named `kt2`. In R, a list is similar
to a python dict AND a python list (since you can look up values by name, or
by index).

The data that 'kt2' holds is named:

```r
    names(kt2)
```

```
## [1] "events"                 "unit_source_statistics"
## [3] "gauge_netcdf_files"     "events_file"           
## [5] "unit_source_file"       "tsunami_events_file"
```
These are like the 'keys' in a python dict -- except we lookup using the `$`
notation (e.g. `kt2$events_file`).

The important variables describing the events are both "data.frames": a table 
of data where the columns can be different data types. In python the package
"pandas" provides something similar.
* `kt2$unit_source_statistics` summarises info on the unit-sources
* `kt2$events` describes the earthquake events, including their slip and occurrence rate, with reference to the unit-sources

Now look at a few rows of the `unit_source_statistics` data.frame.

```r
    head(kt2$unit_source_statistics)
```

```
##      lon_c     lat_c     depth   strike       dip rake slip   length
## 1 186.4512 -14.80783  3.439599 115.6372 10.215011   90    1 54.29835
## 2 186.2351 -15.07410 13.773286 123.6983 20.721931   90    1 48.72223
## 3 186.0175 -15.31206 30.333686 132.7812 29.479209   90    1 39.29113
## 4 186.8934 -15.08891  2.836145 124.5024  7.120315   90    1 59.92791
## 5 186.5826 -15.35606 12.384518 129.4618 16.947532   90    1 48.97650
## 6 186.2618 -15.58282 29.548373 135.6988 26.860916   90    1 40.61590
##      width downdip_number alongstrike_number subfault_number max_depth
## 1 39.20347              1                  1               1  7.477490
## 2 39.00196              2                  1               2 21.553622
## 3 38.72442              3                  1               3 40.000000
## 4 45.51458              1                  2               4  6.280907
## 5 45.35952              2                  2               5 19.781123
## 6 46.10423              3                  2               6 40.000000
##                                                                                                       initial_condition_file
## 1 /g/data1a/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/kermadectonga2/EQ_SOURCE/Unit_source_data/kermadectonga2/kermadectonga2_1_1.tif
## 2 /g/data1a/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/kermadectonga2/EQ_SOURCE/Unit_source_data/kermadectonga2/kermadectonga2_2_1.tif
## 3 /g/data1a/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/kermadectonga2/EQ_SOURCE/Unit_source_data/kermadectonga2/kermadectonga2_3_1.tif
## 4 /g/data1a/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/kermadectonga2/EQ_SOURCE/Unit_source_data/kermadectonga2/kermadectonga2_1_2.tif
## 5 /g/data1a/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/kermadectonga2/EQ_SOURCE/Unit_source_data/kermadectonga2/kermadectonga2_2_2.tif
## 6 /g/data1a/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/kermadectonga2/EQ_SOURCE/Unit_source_data/kermadectonga2/kermadectonga2_3_2.tif
##                                                                                                                                                                                   tide_gauge_file
## 1 /g/data/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/kermadectonga2/TSUNAMI_UNIT_SOURCES/unit_source_tsunami/RUN_20180815163220_kermadectonga2_1_1/RUN_ID100001_20180815_180239.429/Gauges_data_ID100001.nc
## 2 /g/data/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/kermadectonga2/TSUNAMI_UNIT_SOURCES/unit_source_tsunami/RUN_20180815163230_kermadectonga2_2_1/RUN_ID100001_20180815_190101.118/Gauges_data_ID100001.nc
## 3 /g/data/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/kermadectonga2/TSUNAMI_UNIT_SOURCES/unit_source_tsunami/RUN_20180815163245_kermadectonga2_3_1/RUN_ID100001_20180815_200326.773/Gauges_data_ID100001.nc
## 4 /g/data/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/kermadectonga2/TSUNAMI_UNIT_SOURCES/unit_source_tsunami/RUN_20180815163222_kermadectonga2_1_2/RUN_ID100001_20180815_180239.429/Gauges_data_ID100001.nc
## 5 /g/data/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/kermadectonga2/TSUNAMI_UNIT_SOURCES/unit_source_tsunami/RUN_20180815163232_kermadectonga2_2_2/RUN_ID100001_20180815_190931.274/Gauges_data_ID100001.nc
## 6 /g/data/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/kermadectonga2/TSUNAMI_UNIT_SOURCES/unit_source_tsunami/RUN_20180815163249_kermadectonga2_3_2/RUN_ID100001_20180815_200326.797/Gauges_data_ID100001.nc
```
Important info includes:
* `kt2$unit_source_statistics$downdip_number` - a location index, where 1 means "at trench", 2 means "middle row", 3 means "most down-dip"
* `kt2$unit_source_statistics$alongstrike_number` - another location index
* `kt2$unit_source_statistics$subfault_number` - used in `kt2$events` to refer to individual unit-sources.

Next have a look at the events data, which contains 1 row for every event.

```r
    head(kt2$events)
```

```
##   event_index_string       event_slip_string  Mw target_lon target_lat
## 1                 1-                  1.109_ 7.2   186.4512  -14.80783
## 2               1-2-          1.107_0.00229_ 7.2   186.4512  -14.80783
## 3             1-2-4- 1.065_0.03671_0.008215_ 7.2   186.4512  -14.80783
## 4                 1-                  1.109_ 7.2   186.4512  -14.80783
## 5               1-2-          0.9261_0.2044_ 7.2   186.4512  -14.80783
## 6                 1-                  1.109_ 7.2   186.4512  -14.80783
##   peak_slip_downdip_ind peak_slip_alongstrike_ind
## 1                     1                         1
## 2                     1                         1
## 3                     1                         1
## 4                     1                         1
## 5                     1                         1
## 6                     1                         1
##   physical_corner_wavenumber_x physical_corner_wavenumber_y     sourcename
## 1                  0.015456123                  0.023810254 kermadectonga2
## 2                  0.035345239                  0.029880412 kermadectonga2
## 3                  0.020252455                  0.019167928 kermadectonga2
## 4                  0.007980015                  0.016183619 kermadectonga2
## 5                  0.020362803                  0.010794206 kermadectonga2
## 6                  0.012554852                  0.009631231 kermadectonga2
##   uniform_event_row  rate_annual rate_annual_lower_ci rate_annual_upper_ci
## 1                 1 0.0001424545         6.072586e-05         0.0002862334
## 2                 1 0.0001576492         6.720311e-05         0.0003167642
## 3                 1 0.0001788554         7.624292e-05         0.0003593737
## 4                 1 0.0001362682         5.808877e-05         0.0002738034
## 5                 1 0.0002111157         8.999492e-05         0.0004241942
## 6                 1 0.0001307877         5.575251e-05         0.0002627914
##   variable_mu_Mw variable_mu_rate_annual variable_mu_rate_annual_lower_ci
## 1       6.882028            0.0001515197                     6.175740e-05
## 2       6.882843            0.0001528947                     6.231781e-05
## 3       6.894464            0.0001604965                     6.541624e-05
## 4       6.882028            0.0001529224                     6.232911e-05
## 5       6.946053            0.0001780539                     7.257241e-05
## 6       6.882028            0.0001558386                     6.351772e-05
##   variable_mu_rate_annual_upper_ci variable_mu_rate_annual_median
## 1                     0.0003061462                   0.0001458727
## 2                     0.0003089243                   0.0001471964
## 3                     0.0003242839                   0.0001545150
## 4                     0.0003089803                   0.0001472231
## 5                     0.0003597587                   0.0001714181
## 6                     0.0003148725                   0.0001500306
##   variable_mu_rate_annual_16pc variable_mu_rate_annual_84pc
## 1                 0.0000995066                 0.0002169490
## 2                 0.0001004096                 0.0002189177
## 3                 0.0001054019                 0.0002298023
## 4                 0.0001004278                 0.0002189574
## 5                 0.0001169323                 0.0002549413
## 6                 0.0001023429                 0.0002231329
##   variable_mu_weight_with_nonzero_rate weight_with_nonzero_rate
## 1                                    1                        1
## 2                                    1                        1
## 3                                    1                        1
## 4                                    1                        1
## 5                                    1                        1
## 6                                    1                        1
##   rate_annual_16pc rate_annual_84pc rate_annual_median
## 1     8.955377e-05     0.0001995052       0.0001368756
## 2     9.910592e-05     0.0002207851       0.0001514753
## 3     1.124371e-04     0.0002504840       0.0001718510
## 4     8.566480e-05     0.0001908414       0.0001309317
## 5     1.327175e-04     0.0002956640       0.0002028479
## 6     8.221947e-05     0.0001831660       0.0001256658
```
Important columns include:
* `kt2$events$event_index_string` contains the `subfault_number` of each unit-source involved in the event, separated by a `-`.
* `kt2$events$event_slip_string` contains the slip on each of unit-source involved in the event, separated by a `_`.
* `kt2$events$rate_annual` is important, even if you do not care about event frequencies, because **some events are impossible** and
  these have `kt2$events$rate_annual = 0`. The impossible events either have Mw-max being too large, or, they have max-slip being
  too great for their magnitude. See the PTHA18 report for details.

We can get the indices of possible events like so

```r
possible_inds = which(kt2$events$rate_annual > 0)
head(possible_inds)
```

```
## [1] 1 2 3 4 5 6
```

Note you can separate the `event_index_string` and `event_slip_string` into a set of numbers for each event like so:

```r
# Split the strings and cast to numeric. The result will be a list, containing one vector for each event
event_inds_list = lapply( strsplit(kt2$events$event_index_string, '-'), as.numeric)
head(event_inds_list)
```

```
## [[1]]
## [1] 1
## 
## [[2]]
## [1] 1 2
## 
## [[3]]
## [1] 1 2 4
## 
## [[4]]
## [1] 1
## 
## [[5]]
## [1] 1 2
## 
## [[6]]
## [1] 1
```


```r
# Split the slips and cast to numeric. The result will be a list, containing one vector for each event
event_slips_list = lapply( strsplit(kt2$events$event_slip_string, '_'), as.numeric)
head(event_slips_list)
```

```
## [[1]]
## [1] 1.109
## 
## [[2]]
## [1] 1.10700 0.00229
## 
## [[3]]
## [1] 1.065000 0.036710 0.008215
## 
## [[4]]
## [1] 1.109
## 
## [[5]]
## [1] 0.9261 0.2044
## 
## [[6]]
## [1] 1.109
```

## Geting the 3D displacement vector associated with each event.

The events are in a big table named `kt2$events`. The script [get_displacements_for_events.R](get_displacements_for_events.R)
computed their displacements at a chosen location (see `target_pt` therein), and stored it in a matrix named `xyz_displacement_events`.This has 3 columns (x/y/z displacement), and the same number of rows as `kt2$events`. 

Let's print rows 30:50 (corresponding to the events in `kt2$events[30:50,]`)

```r
# First column is easting (m), second is northing (m), third is up/down (m)
xyz_displacement_events[30:50,]
```

```
##               [,1]         [,2]         [,3]
##  [1,] 0.000000e+00 0.0000000000 0.000000e+00
##  [2,] 1.027297e-04 0.0002381513 8.464318e-06
##  [3,] 7.537615e-05 0.0001589102 2.442959e-06
##  [4,] 1.027297e-04 0.0002381513 8.464318e-06
##  [5,] 1.027297e-04 0.0002381513 8.464318e-06
##  [6,] 1.028319e-04 0.0002377538 8.284261e-06
##  [7,] 1.027455e-04 0.0002379095 8.399340e-06
##  [8,] 9.365555e-05 0.0002171153 7.716664e-06
##  [9,] 1.027297e-04 0.0002381513 8.464318e-06
## [10,] 1.015097e-04 0.0002116095 2.719647e-06
## [11,] 1.029704e-04 0.0002368059 7.918607e-06
## [12,] 1.028520e-04 0.0002373814 8.161418e-06
## [13,] 1.027297e-04 0.0002381513 8.464318e-06
## [14,] 1.027297e-04 0.0002381513 8.464318e-06
## [15,] 1.027297e-04 0.0002381513 8.464318e-06
## [16,] 1.027297e-04 0.0002381513 8.464318e-06
## [17,] 0.000000e+00 0.0000000000 0.000000e+00
## [18,] 0.000000e+00 0.0000000000 0.000000e+00
## [19,] 0.000000e+00 0.0000000000 0.000000e+00
## [20,] 0.000000e+00 0.0000000000 0.000000e+00
## [21,] 0.000000e+00 0.0000000000 0.000000e+00
```

Notice how many of the displacements are zero. The reason for this is that in the PTHA18 unit-source construction, we only compute the Okada deformation within a neighbourhood of the unit-source (radius = 20x unit-source depth). This neighbourhood is larger for deep unit-sources, and shallower for near-trench sources. But if earthquakes only include unit-sources far from our `target_pt`, the displacement is zero.

If you would like to use a larger radius for the Okada calculation, to reduce the number of zero displacement events, we can change the variable `okada_distance_factor` in [config.R](config.R) and re-run the unit-source creation code.

Remember that not all of these events are possible according to the PTHA18! This will matter for large events.

# 3) Associate each event and displacement with a tsunami height at, say, Nuku'alofa

