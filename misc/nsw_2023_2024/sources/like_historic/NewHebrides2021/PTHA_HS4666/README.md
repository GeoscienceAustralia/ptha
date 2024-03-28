Scenario identified from earlier work as generally matching New Hebrides 2021 at NSW tide gauges.

Many other scenarios gave a similar match. 

Copied with:
```
cp /media/gareth/Windows7_OS/Users/gareth/Documents/work/Inundation_tsunami/2019_07_australia_wide_model/clean_model_files/ptha18_sce
narios_random/set_range_of_mw_and_centroid_batch3/newhebrides2021-batch3/*4666*.tif .
```

Scenario details
```r
> setwd("/media/gareth/Windows7_OS/Users/gareth/Documents/work/Inundation_tsunami/2019_07_australia_wide_model/clean_model_files/ptha18_scenarios_random/set_range_of_mw_and_centroid_batch3")

> load("find_scenarios_near_historical_events_R_image_C.RData")

# We want 4666 -- index 51
> random_scenarios[["newhebrides2021-batch3"]][['heterogeneous_slip']]$desired_event_rows
 [1] 5688 6336 4652 5683 5682 5730 4593 6330 5696 5703 6414 6280 4706 5728 6383
[16] 4736 6244 6300 5698 5717 6297 6324 5726 4738 5712 6344 4711 4708 6304 4665
[31] 5705 5671 4701 4733 4688 6384 5674 4693 5715 6325 6237 5690 4696 4714 4699
[46] 6242 4686 4679 4683 4623 4666 6261 5741 5706 4734 4635
> which(random_scenarios[["newhebrides2021-batch3"]][['heterogeneous_slip']]$desired_event_rows == 4666)
[1] 51

> random_scenarios[["newhebrides2021-batch3"]][['heterogeneous_slip']]$events[51,]
   event_index_string   event_slip_string  Mw target_lon target_lat
51          12-14-16- 0.7422_2.668_1.504_ 7.6   171.4139  -22.39378
   peak_slip_downdip_ind peak_slip_alongstrike_ind physical_corner_wavenumber_x
51                     2                         7                  0.004993983
   physical_corner_wavenumber_y   sourcename uniform_event_row rate_annual
51                  0.008211274 newhebrides2               312 2.21587e-05
   rate_annual_lower_ci rate_annual_upper_ci variable_mu_Mw
51         2.702381e-05         4.255709e-05        7.72709
   variable_mu_rate_annual variable_mu_rate_annual_lower_ci
51            2.176725e-05                     2.611479e-05
   variable_mu_rate_annual_upper_ci variable_mu_rate_annual_median
51                     4.620979e-05                   1.904445e-05
   variable_mu_rate_annual_16pc variable_mu_rate_annual_84pc
51                 1.325878e-05                 2.920542e-05
   variable_mu_weight_with_nonzero_rate weight_with_nonzero_rate
51                                    1                        1
   rate_annual_16pc rate_annual_84pc rate_annual_median
51     1.367234e-05     2.958419e-05       1.918197e-05

```
