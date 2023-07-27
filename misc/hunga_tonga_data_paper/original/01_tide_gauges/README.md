# Tide gauge observations in Australia during the January 2022 tsunami

This folder contains observations at a number of tide-gauges, from a number of data providers. The formats and time-zones vary. 

* [2022_Tsunami_TAS](2022_Tsunami_TAS) contains data for sites in Tasmania, collated by Karen Palmer during her PhD. See metadata in the header of each file. 
    * The data is recorded every 6 min using pressure sensors in a stilling well, and later corrected for atmospheric pressure using nearby pressure sensors. 
    * The stilling well pressure is an instantaneous measurement, although the stilling well will cause some mechanical smoothing.
* [AAD_2022_tidegauge](AAD_2022_tidegauge) contains tide gauge data for Davis, Mawson, and Casey, provided by the Australian Antarctic Division.
    * The tide gauges measure the sum of the atmospheric and sea level pressure.
    * The records have a 10 min spacing. They represent a 10 min period of continuously averaged 1s data. This will tend to damp the tsunami signal.
    * The folder also contains nearby MSLP measurements. They are interpolated to the tide gauges and subtracted from the tide gauge pressure to estimate the pressure due to sea level variations alone. This assumes a pressure change of 1 dbar corresponds to a sea level change of 1 m.
* [BOM_tidegaugedata](BOM_tidegaugedata) contains data for a range of sites in Australia, collated by the Bureau of Meteorology as part of collaborations with various Port Authorities (BOMPorts). 
    * For metadata see the file [BOM_and_MHL_tide_gauge_metadata_Table.csv](BOM_and_MHL_tide_gauge_metadata_Table.csv). 
    * The instruments and sampling regimes are variable, and we do not have details for these stations.
* [DES_QGHL_data](DES_QGHL_data) contains data for sites in Queensland. 
    * See metadata in the header of each file. 
    * These gauges represent a 4-minute average of the waterlevel (this is done internally in the instrument, before storing).
* [ioc_sealevelmonitoring](ioc_sealevelmonitoring) contains data for Australian stations that is available from the [IOC sea level monitoring website](https://www.ioc-sealevelmonitoring.org/). 
    * For metadata see the latter website, as well as the file [ioc_sealevelmonitoring/ioc_australian_station_list.csv](ioc_sealevelmonitoring/ioc_australian_station_list.csv). 
    * These stations are maintained by the Bureau of Meteorology.
* [Macquarie_Island](Macquarie_Island) contains data and metadata for a tide gauge on Macquarie Island, provided by the Australian Antarctic Division. 
    * The records have 3 min spacing, and represent an average of 1 s measurements over that period.
* [MHL_data](MHL_data) contains data for a range of sites in New South Wales (and also Lord Howe Island), collected by Manly Hydraulics Laboratory. 
    * For metadata see [the pdf file inside that folder](MHL_data/Data conditions and limitations.pdf) and also the file [BOM_and_MHL_tide_gauge_metadata_Table.csv](BOM_and_MHL_tide_gauge_metadata_Table.csv).
    * MHLs one minute tide gauge data represents the average of 60 1s samples, taken 30s either side of the reported time.
* [NSW_Port_Authority_DO_NOT_DISTRIBUTE](NSW_Port_Authority_DO_NOT_DISTRIBUTE) is an empty folder because it contained data that we do not have permission to distribute, collected by the Port Authority of NSW. 
    * However, we do have permission to distribute de-tided versions of the data. 
    * In practice the restricted data was stored in this folder for processing, and then removed prior to distribution. 
    * The tide gauges employ a range of sampling methods, and we do not have details for all stations. 
        * At Fort Denison, Sydney, an acoustic pulse measures water levels in a stilling well every second. Measurements are derived by averageing 60 of these over 1 minute. 
        * At the Cruise Wharf, Eden, there is no stilling well and a Vega radar gauge measures an average water level every 5 seconds, with 12 of these averaged to produce the 1 minute data. In addition the instrument uses an in-built algorithm to smooth out large fluctiotions between readings, which can attenuate large water level shifts with short periods.
