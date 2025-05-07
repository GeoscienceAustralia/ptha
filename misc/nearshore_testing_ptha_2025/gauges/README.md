# Uniform interface to a range of detided tide gauge data
---------------------------------------------------------

The script `gauge_data_links.R` provides a uniform interface to detided tide-gauge data that we have obtained from many sources over the years. It contains a list of tide gauges with:
* Location information
* A list of files containing tide gauge observations for different tsunami events. These were created by the author using datasets that were obtained as described in the `./DATA/` folder near each dataset (download as explained below).
* A function which reads this data and converts it to a uniform format. This ensures a uniform time zone (UTC) and consistent column names, which the underlying data does not have.

To use the script you have to download and extract the associated data files here: https://thredds.nci.org.au/thredds/fileServer/fj6/PTHA/Nearshore_testing_2025/DATA.tar.bz2
* They can be extracted with (e.g.) `tar -jxf DATA.tar.bz2`
* You should end up with a folder `./DATA/` containing all the files that `gauge_data_links.R` points to.

Note that code inside `../swals` uses a different file path to refer to the `gauge_data_links.R` script, since on the machines used to run the analysis, the tide gauge data is stored in a central location (not here). 


## Usage
--------

Start R and source the `gauge_data_links.R` script.
```r
source('gauge_data_links.R') 
```

Each tsunami is associated with a string giving the day it occurred in UTC time.
  * For example, the 2021 Kermadec tsunami occurred on March 4 2021, and uses `'2021-03-04'`

Get the data for that event using `get_data_for_event`
```r
gauge_data_kermadec2021 = get_data_for_event('2021-03-04')
```

The result is a named list with one entry for each gauge with data for the
event (26 gauges here). 
```r
length(gauge_data_kermadec2021)
# [1] 26

names(gauge_data_kermadec2021)
#  [1] "BallinaBreakwall_1min_DPIE"             
#  [2] "Bermagui_1min_DPIE"                     
#  [3] "BrunswickHeads_1min_DPIE"               
#  [4] "CoffsHarbourInnerPumpoutJetty_1min_DPIE"
#  [5] "CrookhavenHeads_1min_DPIE"              
#  [6] "CrowdyHeadFishermansWharf_1min_DPIE"    
#  [7] "Forster_1min_DPIE"                      
#  [8] "LordHoweIsland_1min_DPIE"               
#  [9] "PortMacquarie_1min_DPIE"                
# [10] "ShoalBay_1min_DPIE"                     
# [11] "TweedEntranceSouth_1min_DPIE"           
# [12] "Yamba_1min_DPIE"                        
# [13] "Eden_1min_DPIE"                         
# [14] "TwofoldBay_1min_BOM"                    
# [15] "GoldCoastSandBypass"                    
# [16] "BatemansBay_PrincessJetty_1min_DPIE"    
# [17] "Ulladulla_1min_DPIE"                    
# [18] "JervisBay_1min_DPIE"                    
# [19] "Bundeena_1min_DPIE"                     
# [20] "Sydney_MiddleHarbour_1min_DPIE"         
# [21] "Sydney_FortDenison_1min_PA_b"           
# [22] "Sydney_Botany_Bay_Pilot_Jetty_1min_PA"  
# [23] "Hawkesbury_Patonga_1min_DPIE"           
# [24] "Newcastle_east_1min_PA"                 
# [25] "PortKembla_1min_BOM"                    
# [26] "PortKembla_GrainTerminal_1min_PA"       
```

For each gauge we can extract the key information like this (using `'Eden_1min_DPIE'` as an example):
```r
gn = 'Eden_1min_DPIE' # Gauge name

gauge_data_kermadec2021[[gn]]$coord # Coordinate (lon, lat)
# [1] 149.90829 -37.07124


site_data = gauge_data_kermadec2021[[gn]]$obs # this is a data.frame 

options(digits=12) # Reduce rounding of printed output below
head(site_data) # Print the first few rows
#                  time       juliant height             resid
# 1 2021-01-31 14:00:00 18658.5833333  1.053 -0.01150296973029
# 2 2021-01-31 14:01:00 18658.5840278  1.059 -0.00293755044633
# 3 2021-01-31 14:02:00 18658.5847222  1.059 -0.00036937101384
# 4 2021-01-31 14:03:00 18658.5854167  1.053 -0.00379850327846
# 5 2021-01-31 14:04:00 18658.5861111  1.048 -0.00622501968654
# 6 2021-01-31 14:05:00 18658.5868056  1.033 -0.01864899320414
```
The data contains the observation `time` (UTC), `juliant` (time in days
since the start of 1970 UTC), `height` (the tide gauge value with an arbitrary
datum) and `resid` (the detided record).

The script ensures that all sites are in this format. If you just want to work with
the csv files, they are in the `./DATA` archive, but don't have a consistent timezone
or columns or coordinates. However, you can read `gauge_data_links.R` to figure out this
information.
