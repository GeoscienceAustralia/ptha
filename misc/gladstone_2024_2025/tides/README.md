# Tidal Corrections

The region between Yeppoon and Burnet heads has a varying HAT level. This folder contains the data and scripts to correct the tides for this region.

## Data
- Tidal model predictions stored in /g/data/w85/tsunami/MODELS/tidal_prediction/TPXO9
- Tide gauge records of HAT stored in `x_y_z.csv` compiled from online sources:
    - QLD storm tide website, which comes from the tide tables. https://www.qld.gov.au/environment/coasts-waterways/beach/storm/storm-sites.

## Scripts
Originally I startd making `interpolate_hat.R` (deleted) as we weren't going to use TPXO9 for licencing reasons, but we have since been given permission to use it. The script `adjust_tpxo9_for_gauges.R` is currently being used to adjust the TPXO9 predictions to the tide gauge records region.

## Outputs
The resulting geotif file will be read by the model to adjust the tides for the region. The model will add the geotif to the elevation data, so the results in `tpxo_adjusted_for_gauges.tif` are negated to add depth to the sea level.

## Tide Gauge Archives
Historical tide gauge data were downloaded for selected sites from https://www.msq.qld.gov.au/Tides/Open-data and placed into the [tide_gauge_archives](tide_gague_archives) directory. The script [tide_gauge_archives/hat_exceedence.R](tide_gauge_archives/hat_exceedance.R) plots this data as violin plot grouped by year for each site, the timeseries and the fraction of observations that exceed a given threshold (specified in the script), i.e. HAT (m AHD).

## Tests
1. Run `test_adjustment.R` to ensure that the geotif file created does indeed match the tide gagues.
Pass
2. Compare model runs with a static tide and a HAT adjusted tide.
This was done for three locations: Lady elliot Island (HAT=1.41 m AHD), Rosslyn Bay (HAT = 2.85 m AHD) and South Trees (HAT = 2.42 m AHD).
