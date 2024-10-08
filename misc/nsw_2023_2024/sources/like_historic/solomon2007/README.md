This source model is from PTHA18
* Solomon 2 source, VAUS slip model, row index 11244

My unpublished studies using an Australia wide model (setup as per the Global Dissipation Models paper) suggest it well simulates tide gauge observations in SE Australia.

Extracted from here: `"/media/gareth/Windows7_OS/Users/gareth/Documents/work/Inundation_tsunami/2019_07_australia_wide_model/clean_model_files/ptha18_scenarios_random/set_range_of_mw_and_centroid_batch2"`

Some summary statistiscs are below.
```r

> # Running in this folder
> getwd()
[1] "/media/gareth/Windows7_OS/Users/gareth/Documents/work/Inundation_tsunami/2019_07_australia_wide_model/clean_model_files/ptha18_scenarios_random/set_range_of_mw_and_centroid_batch2"

> # Get the scenario summary stats relevant to this event
> load('find_scenarios_near_historical_events_R_image_C.RData')

# We are looking for 11244 (I know that from the plots). 
# Here we see it corresponds to index 26 in the scenarios
> random_scenarios[["solomon2007-batch2"]][["variable_area_uniform_slip"]]$desired_event_rows
 [1] 10091 10029 10657 11174 10606 10640 10664 10678 10059 10684 10588 10076
[13] 11200 10623 11264 10676 10047 10611 11213  9991  9982 10654 10037 10645
[25] 10625 11244 11225 11261 10054

# Here are the summary statistics
> random_scenarios[["solomon2007-batch2"]][["variable_area_uniform_slip"]]$events[26,]
   event_index_string        event_slip_string  Mw target_lon target_lat
26       39-40-41-42- 8.739_8.739_8.739_8.739_ 8.2   155.3399  -6.939554
   peak_slip_downdip_ind peak_slip_alongstrike_ind physical_corner_wavenumber_x
26                     1                        21                  0.004478087
   physical_corner_wavenumber_y sourcename uniform_event_row rate_annual
26                   0.01089278   solomon2               750 3.37831e-05
   rate_annual_lower_ci rate_annual_upper_ci variable_mu_Mw
26         5.126155e-06         5.222016e-05       8.291904
   variable_mu_rate_annual variable_mu_rate_annual_lower_ci
26            4.495603e-05                     8.554527e-06
   variable_mu_rate_annual_upper_ci variable_mu_rate_annual_median
26                     7.340527e-05                   4.072247e-05
   variable_mu_rate_annual_16pc variable_mu_rate_annual_84pc
26                 4.000076e-05                 5.439946e-05
   variable_mu_weight_with_nonzero_rate weight_with_nonzero_rate
26                            0.9799688                0.9789364
   rate_annual_16pc rate_annual_84pc rate_annual_median
26     3.485496e-05     4.019378e-05       3.057462e-05
```
