Scenario identified from earlier work as matching NSW Sumatra 2004.

Copied with
```
cp /media/gareth/Windows7_OS/Users/gareth/Documents/work/Inundation_tsunami/2019_07_australia_wide_model/clean_model_files/ptha18_scenarios_random/set_range_of_mw_and_centroid_batch2/sumatra2004-batch2/*107476* .
```

Scenario details
```r
> setwd("/media/gareth/Windows7_OS/Users/gareth/Documents/work/Inundation_tsunami/2019_07_australia_wide_model/clean_model_files/ptha18_scenarios_random/set_range_of_mw_and_centroid_batch2")

> load("find_scenarios_near_historical_events_R_image_C.RData")

# We want 107476 -- this is in index 4
> random_scenarios[["sumatra2004-batch2"]]$variable_area_uniform_slip$desired_event_rows
 [1] 104599 104412 107575 107476 104284 106059 105982 105877 104394 107570
[11] 104674 104477 106119 104417 104691 107364 105876 104548 105850 107665
[21] 104294 106121 104628 104623 107023 104395 105886 104406 106070

# Data
> random_scenarios[["sumatra2004-batch2"]]$variable_area_uniform_slip$events[4,]
                                                                                                                                                event_index_string
4 277-278-281-282-285-286-289-290-293-294-297-298-301-302-305-306-309-310-313-314-317-318-321-322-325-326-329-330-333-334-337-338-341-342-345-346-349-350-353-354-
                                                                                                                                                                                                                                 event_slip_string
4 28.76_28.76_28.76_28.76_28.76_28.76_28.76_28.76_28.76_28.76_28.76_28.76_28.76_28.76_28.76_28.76_28.76_28.76_28.76_28.76_28.76_28.76_28.76_28.76_28.76_28.76_28.76_28.76_28.76_28.76_28.76_28.76_28.76_28.76_28.76_28.76_28.76_28.76_28.76_28.76_
   Mw target_lon target_lat peak_slip_downdip_ind peak_slip_alongstrike_ind
4 9.3   94.31102   5.261888                     1                        74
  physical_corner_wavenumber_x physical_corner_wavenumber_y sourcename
4                  0.001392028                  0.003432553     sunda2
  uniform_event_row  rate_annual rate_annual_lower_ci rate_annual_upper_ci
4              7166 7.115709e-07                    0           1.5575e-06
  variable_mu_Mw variable_mu_rate_annual variable_mu_rate_annual_lower_ci
4       9.004051            5.384678e-07                                0
  variable_mu_rate_annual_upper_ci variable_mu_rate_annual_median
4                     1.212523e-06                   7.481909e-07
  variable_mu_rate_annual_16pc variable_mu_rate_annual_84pc
4                            0                 9.339598e-07
  variable_mu_weight_with_nonzero_rate weight_with_nonzero_rate
4                            0.7005843                0.7248933
  rate_annual_16pc rate_annual_84pc rate_annual_median
4                0     1.170948e-06        9.20221e-07
```
