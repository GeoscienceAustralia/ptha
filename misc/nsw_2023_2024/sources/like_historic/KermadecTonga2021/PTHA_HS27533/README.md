Scenario identified from earlier work as generally matching Kermadec 2021 at NSW tide gauges.

Many other scenarios gave a similar match. 

Copied with:
```
cp /media/gareth/Windows7_OS/Users/gareth/Documents/work/Inundation_tsunami/2019_07_australia_wide_model/clean_model_files/ptha18_scenarios_random/set_range_of_mw_and_centroid_batch3/kermadec2021-batch3/*27533* .
```

Scenario details
```r
> setwd("/media/gareth/Windows7_OS/Users/gareth/Documents/work/Inundation_tsunami/2019_07_australia_wide_model/clean_model_files/ptha18_scenarios_random/set_range_of_mw_and_centroid_batch3")

> load("find_scenarios_near_historical_events_R_image_C.RData")

# We want 27533
> which(random_scenarios[["kermadec2021-batch3"]][['heterogeneous_slip']]$desired_event_rows == 27533)
[1] 45

>  random_scenarios[["kermadec2021-batch3"]][['heterogeneous_slip']]$events[45,]
             event_index_string                           event_slip_string  Mw
45 103-105-106-107-108-109-111- 2.49_3.471_8.235_0.06115_3.065_5.758_1.441_ 8.1
   target_lon target_lat peak_slip_downdip_ind peak_slip_alongstrike_ind
45  -176.5247  -28.79711                     1                        36
   physical_corner_wavenumber_x physical_corner_wavenumber_y     sourcename
45                  0.003678311                  0.007859026 kermadectonga2
   uniform_event_row  rate_annual rate_annual_lower_ci rate_annual_upper_ci
45              1836 4.135299e-06         2.399478e-06         6.975115e-06
   variable_mu_Mw variable_mu_rate_annual variable_mu_rate_annual_lower_ci
45       8.020481            4.061399e-06                     2.638988e-06
   variable_mu_rate_annual_upper_ci variable_mu_rate_annual_median
45                      6.75034e-06                   3.654398e-06
   variable_mu_rate_annual_16pc variable_mu_rate_annual_84pc
45                 2.751663e-06                 5.059869e-06
   variable_mu_weight_with_nonzero_rate weight_with_nonzero_rate
45                                    1                        1
   rate_annual_16pc rate_annual_84pc rate_annual_median
45     2.810033e-06     5.237954e-06       3.846768e-06

