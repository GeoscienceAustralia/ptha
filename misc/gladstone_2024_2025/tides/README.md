# Tidal Corrections

The region between Yeppoon and Burnet heads has a varying HAT level. This folder contains the data and scripts to correct the tides for this region.

## Data
- Tidal model predictions stored in /g/data/w85/tsunami/MODELS/tidal_prediction/TPXO9
    - Note we were given written permission by the TPXO9 authors to use it for this study
- Tide gauge records of HAT stored in `x_y_z.csv` compiled from online sources:
    - QLD storm tide website, which comes from the tide tables. https://www.qld.gov.au/environment/coasts-waterways/beach/storm/storm-sites.

## Scripts
The QLD regional raster extracted from TPXO9 of the highest astronomical tides is treated in [tpxo9_adjusted_for_gauges.R](tpxo9_adjusted_for_gauges.R).
1. First onshore areas where TPXO9 is empty are extrapolated onto by filling westwards from the westernmost TPXO9 value. This avoids an abrupt transition at the coastline in the tidal adjustment. 
2. Then a background sea level is selected from the mean of the cells on the eastern side of the TPXO9 raster. Then the values around the edge of raster are linearly transitioned to the background sea level over a 1 degree transition zone (within the raster).
3. The raster is extended to match the SWALS model extent, filing in with the background sea level.

There's also a function to smoothly bump the sea level up or down locally. Although this is not used. Evidence from tide gauge records could suggest local variations from TPXO9 in which case this could be used.

The resulting geotif file will be read by the model to adjust the tides for the region. The model will add the geotif to the elevation data, so the results in `tpxo_adjusted_for_gauges.tif` are negated to add depth to the sea level.

## Tide Gauge Archives
Historical tide gauge data were downloaded for selected sites from https://www.msq.qld.gov.au/Tides/Open-data and placed into the [tide_gauge_archives](tide_gauge_archives) directory. The script [tide_gauge_archives/hat_exceedance.R](tide_gauge_archives/hat_exceedance.R) plots this data as violin plot grouped by year for each site, the timeseries and the fraction of observations that exceed a given threshold (specified in the script), i.e. HAT (m AHD).

## Tests
1. Run `test_adjustment.R` to ensure that the geotif file created does indeed match the tide gagues.
Pass
2. Compare model runs with a static tide and a HAT adjusted tide.
This was done for three locations: Lady elliot Island (HAT=1.41 m AHD), Rosslyn Bay (HAT = 2.85 m AHD) and South Trees (HAT = 2.42 m AHD).
