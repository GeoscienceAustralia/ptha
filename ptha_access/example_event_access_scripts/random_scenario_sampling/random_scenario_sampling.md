# Randomly sample PTHA18 scenarios on a source-zone
-------------------------------------------------

The PTHA18 often includes thousands or tens-of-thousands of scenarios on a
source-zone. For some applications it is impractical to work with all
scenarios, but may be practical to work with a random sample of scenarios.

This tutorial shows how to randomly sample scenarios from a given source-zone in
a manner that respects the PTHA18 scenario conditional probabilities and earthquake
magnitude-frequency modelling. 

## Get the source-zone event data
---------------------------------

Here we choose to work with heterogeneous-slip scenarios from the
`kermadectonga2` source-zone. 

Firstly we source the scripts to work with the PTHA18 data.

```r
# Get the scripts to access the PTHA18
ptha18 = new.env()
source('../../get_PTHA_results.R', local=ptha18, chdir=TRUE)
```

Next we read all scenarios on the kermadectonga2 source-zone.
Here `kt2_scenarios` is a list containing the event table, 
the unit-source-statistics table, and the names of various files
holding data for the source-zone

```r
# Read all heterogeneous-slip scenario metadata (slip_type='stochastic' in
# PTHA18)
kt2_scenarios = ptha18$get_source_zone_events_data(
    'kermadectonga2', 
    slip_type='stochastic')
names(kt2_scenarios)
```

```
## [1] "events"                 "unit_source_statistics" "gauge_netcdf_files"    
## [4] "events_file"            "unit_source_file"       "tsunami_events_file"
```

There are over 44-thousand scenarios in the database, with one
row per scenario in the event table

```r
# How many rows in the scenario events table?
nrow(kt2_scenarios$events)
```

```
## [1] 44685
```

```r
# What are the column-names of the event table?
names(kt2_scenarios$events)
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

For our purposes, one important variable is `Mw` (scenario moment magnitude).
The moment magnitudes are binned and cover the range 7.2, 7.3, ..., 9.8
although the higher magnitudes may be impossible (in which case all scenarios
are assigned a rate of 0).

```r
# Print all unique Mw values -- beware not all of these will be "possible" according
# to PTHA18, as some will have zero probability of occurring.
unique(kt2_scenarios$events$Mw)
```

```
##  [1] 7.2 7.3 7.4 7.5 7.6 7.7 7.8 7.9 8.0 8.1 8.2 8.3 8.4 8.5 8.6 8.7 8.8 8.9 9.0
## [20] 9.1 9.2 9.3 9.4 9.5 9.6 9.7 9.8
```
In this tutorial another important variable is the `rate_annual`, which varies for
each scenario. This is equal to the scenario conditional probability (conditional on
the occurrance of an earthquake with the same magnitude) multiplied by the
logic-tree-mean-rate of scenarios with that magnitude (events/year), according
to the PTHA18. 

Recall that the scenario conditional probability varies between scenarios, and
is used to account for spatial variations in tectonic convergence, to limit the
scenraio peak-slip, and to adjust for bias in the earthquake-source models
(although this is more prominent for variable-area-uniform-slip scenarios than
for heterogeneous-slip scenarios).  See 
[this paper](https://doi.org/10.1007/s00024-019-02299-w) for further information.

## Random scenario sampling, stratified by magnitude
----------------------------------------------------

Our simplest random scenario sampling algorithm proceeds as follows
* Group the scenarios by magnitude
* For each magnitude, sample a given number of scenarios with replacement, with the chance of sampling each scenario proportional to its conditional probability.

This can be implemented as follows (herein we select 12 scenarios for each magnitude)

```r
# Convenient shorthand for the magnitudes and rates in the event table
event_Mw = kt2_scenarios$events$Mw
event_rates = kt2_scenarios$events$rate_annual

# Make the random scenarios
random_scenarios_simple = ptha18$randomly_sample_scenarios_by_Mw_and_rate(
    event_rates=event_rates,
    event_Mw=event_Mw,
    samples_per_Mw=function(Mw){ 12 },
    mw_limits=c(7.15, 9.85))
```

The result is a `data.frame` containing the indices of the random scenarios `inds`,
their magnitudes, `mw`, as well as information on the scenario rates that will be discussed
further below.


```r
# Look at the first few rows
head(random_scenarios_simple)
```

```
##   inds  mw rate_with_this_mw importance_sampling_scenario_rates
## 1  651 7.2        0.05704921                        0.004754101
## 2  504 7.2        0.05704921                        0.004754101
## 3  413 7.2        0.05704921                        0.004754101
## 4  972 7.2        0.05704921                        0.004754101
## 5 1464 7.2        0.05704921                        0.004754101
## 6  235 7.2        0.05704921                        0.004754101
##   importance_sampling_scenario_rates_self_normalised
## 1                                        0.004754101
## 2                                        0.004754101
## 3                                        0.004754101
## 4                                        0.004754101
## 5                                        0.004754101
## 6                                        0.004754101
##   importance_sampling_scenario_weights
## 1                           0.08333333
## 2                           0.08333333
## 3                           0.08333333
## 4                           0.08333333
## 5                           0.08333333
## 6                           0.08333333
##   importance_sampling_scenario_weights_self_normalised
## 1                                           0.08333333
## 2                                           0.08333333
## 3                                           0.08333333
## 4                                           0.08333333
## 5                                           0.08333333
## 6                                           0.08333333
```
Here `inds` corresponds to indices in the `event_Mw` and `event_rates`
variables. Because these were created directly from the event table, the event
indices also correspond to rows in `kt2_scenarios$events`.


In PTHA18 some earthquake magnitudes are impossible. In this case the scenario index will
take an `NA` value, as will various other variables. We see this at the end of the current
table, for magnitudes `9.7` and `9.8`.


```r
# Look at the last few rows
tail(random_scenarios_simple)
```

```
##      inds  mw rate_with_this_mw importance_sampling_scenario_rates
## 297 44105 9.6      5.323646e-05                       4.436371e-06
## 298 44268 9.6      5.323646e-05                       4.436371e-06
## 299 44100 9.6      5.323646e-05                       4.436371e-06
## 300 44138 9.6      5.323646e-05                       4.436371e-06
## 301    NA 9.7      0.000000e+00                                 NA
## 302    NA 9.8      0.000000e+00                                 NA
##     importance_sampling_scenario_rates_self_normalised
## 297                                       4.436371e-06
## 298                                       4.436371e-06
## 299                                       4.436371e-06
## 300                                       4.436371e-06
## 301                                                 NA
## 302                                                 NA
##     importance_sampling_scenario_weights
## 297                           0.08333333
## 298                           0.08333333
## 299                           0.08333333
## 300                           0.08333333
## 301                                   NA
## 302                                   NA
##     importance_sampling_scenario_weights_self_normalised
## 297                                           0.08333333
## 298                                           0.08333333
## 299                                           0.08333333
## 300                                           0.08333333
## 301                                                   NA
## 302                                                   NA
```

Aside from the impossible magnitudes, we can confirm that we have 12 scenarios per magnitude, as requested.

```r
table(random_scenarios_simple$mw)
```

```
## 
## 7.2 7.3 7.4 7.5 7.6 7.7 7.8 7.9   8 8.1 8.2 8.3 8.4 8.5 8.6 8.7 8.8 8.9   9 9.1 
##  12  12  12  12  12  12  12  12  12  12  12  12  12  12  12  12  12  12  12  12 
## 9.2 9.3 9.4 9.5 9.6 9.7 9.8 
##  12  12  12  12  12   1   1
```
