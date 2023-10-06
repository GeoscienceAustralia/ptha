# Mean sea-level pressure data from the Bureau of Meterology and QLD Department of Environment and Science

This folder contains:

* DES pressure station locations reported in [DES_station_locations_all.csv](DES_station_locations_all.csv)
    * Some site names differ slightly from those used in the associated time series files, but it's clear how they match up.

* DES pressure data for each station in [DES_QGHL_pressure_data](DES_QGHL_pressure_data).

* Mean sea level pressure data for each station in [BOM_mslpdata](BOM_mslpdata)
  * Each file contains 3 columns with the `time (timezone = UTC)`, `Mean sea level pressure (hectopascals)`, `Quality Flag (always Y for this data)`

* Metadata for Barometric pressure stations in [HD01D_StnDet.csv](HD01D_StnDet.csv). The columns are:
    1.  Record identifier - st
    2.  Bureau of Meteorology Station Number
    3.  Rainfall district code
    4.  Station Name
    5.  Month/Year station opened (MM/YYYY)
    6.  Month/Year station closed if applicable (MM/YYYY)
    7.  Latitude to 4 decimal places, in decimal degrees
    8.  Longitude to 4 decimal places, in decimal degrees
    9.  Method by which latitude/longitude was derived
    10. State
    11. Height of station above mean sea level in metres
    12. Height of barometer above mean sea level in metres
    13. WMO (World Meteorological Organisation) Index Number
    14. First year of data supplied in data file
    15. Last year of data supplied in data file
    16. Percentage complete between first and last records
    17. Percentage of values with quality flag 'Y' (quality controlled and acceptable)
    18. Percentage of values with quality flag 'N' (not quality controlled)
    19. Percentage of values with quality flag 'W' (quality controlled and considered wrong)
    20. Percentage of values with quality flag 'S' (quality controlled and considered suspect)
    21. Percentage of values with quality flag 'I' (quality controlled and inconsistent with other known information)
    22. `#` symbol, end of record indicator

