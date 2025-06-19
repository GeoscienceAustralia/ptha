# Sources for Solomon 2007 event

The four earthquake unit sources in PTHA18 closest to the global CMT reported centroid were fetched using `fetch_ptha_sources.R`. These are in the `ptha_unit_source` folder. The script then scales them to by magnitude 8.1 and saves them in `ptha_sources_scaled` with their scaling factor. The scaled `solomon2_1_19` event works best for the tide gauges near Gladstone, which is what is used. The unit source has unit source properties:
```
lon_c = 156.483168145658
lat_c = -8.06934791682661
depth = 9.3470115706021
strike = 318.506625534382
dip = 25.2188588328763
rake = 90
slip = 1
length = 54.1742348317985
width = 43.7762585418987
downdip_number = 1
alongstrike_number = 19
subfault_number = 37
max_depth = 20.0640617114261
initial_condition_file = "/g/data1a/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/solomon2/EQ_SOURCE/Unit_source_data/solomon2/solomon2_1_19.tif"
tide_gauge_file = "/g/data/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/solomon2/TSUNAMI_UNIT_SOURCES/unit_source_tsunami/RUN_20180815145226_solomon2_1_19/RUN_ID100001_20180815_145823.782/Gauges_data_ID100001.nc"
distance = 0.313898672450932
```
which is then scaled to a Magnitude 8.1.

A couple of other sources were tried: a) Solomon 2 VAUS row index 11244 from PTHA18, b) Wei et al., 2015, Pure and Applied Geophysics. Sources https://sift.pmel.noaa.gov/ComMIT/compressed//nv_011_b_ha.nc with a scaling factor of alpha = 7.924.
(and `linCo_source001h.nc`) from NOAA had uneven grid spacing, making them difficult to input to SWALS.
