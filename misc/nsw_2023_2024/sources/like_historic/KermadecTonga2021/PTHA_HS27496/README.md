Scenario identified from earlier work as generally matching Kermadec 2021 at NSW tide gauges.

Many other scenarios gave a similar match. 

Copied with:
```
cp /media/gareth/Windows7_OS/Users/gareth/Documents/work/Inundation_tsunami/2019_07_australia_wide_model/clean_model_files/ptha18_scenarios_random/set_range_of_mw_and_centroid_batch3/kermadec2021-batch3/*27496* .
```

Scenario details
```r
> setwd("/media/gareth/Windows7_OS/Users/gareth/Documents/work/Inundation_tsunami/2019_07_australia_wide_model/clean_model_files/ptha18_scenarios_random/set_range_of_mw_and_centroid_batch3")

> load("find_scenarios_near_historical_events_R_image_C.RData")

# We want 27496
> which(random_scenarios[["kermadec2021-batch3"]][['heterogeneous_slip']]$desired_event_rows == 27496)
[1] 21

> random_scenarios[["kermadec2021-batch3"]][['heterogeneous_slip']]$events[21,]
                                     event_index_string
21 95-97-98-99-100-101-102-103-104-105-106-107-108-110-
                                                                     event_slip_string
21 2.192_0.9831_3.591_3.001_2.053_1.794_3.047_1.745_2_1.772_0.3308_1.367_1.046_0.8445_
    Mw target_lon target_lat peak_slip_downdip_ind peak_slip_alongstrike_ind
21 8.1  -176.4096  -28.38701                     2                        33
   physical_corner_wavenumber_x physical_corner_wavenumber_y     sourcename
21                  0.004537711                  0.005947553 kermadectonga2
   uniform_event_row  rate_annual rate_annual_lower_ci rate_annual_upper_ci
21              1834 6.216658e-06         3.607172e-06          1.04858e-05
   variable_mu_Mw variable_mu_rate_annual variable_mu_rate_annual_lower_ci
21       8.071464             5.18944e-06                     3.371958e-06
   variable_mu_rate_annual_upper_ci variable_mu_rate_annual_median
21                     8.625226e-06                   4.669395e-06
   variable_mu_rate_annual_16pc variable_mu_rate_annual_84pc
21                 3.515928e-06                 6.465231e-06
   variable_mu_weight_with_nonzero_rate weight_with_nonzero_rate
21                                    1                        1
   rate_annual_16pc rate_annual_84pc rate_annual_median
21     4.224365e-06     7.874296e-06       5.782904e-06
