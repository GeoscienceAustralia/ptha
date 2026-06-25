# Outputs from Geraldton to Coral Bay, WA tsunami inundation modelling project, 2025-2028

This folder contains outputs for the "first half" of the modelling in Phase 4
of the WA tsunami inundation modelling project. Results cover the region from
Geraldton to Coral Bay inclusive (although results at Geraldton are mainly provided
for comparison with Phase 1-3).

Later, the second half of the modelling will cover the region from Coral Bay to
Onslow, and we will also develop a technical report and provide a
combined-model data package.

## Structure of the model output folders

All files are loaded in a QGIS session contained in the folder
`./QGIS_SESSION`. The same folder contains a backup copy of the QGIS session.

If the QGIS session doesn't work for you, or you'd rather use them in another
way, information on the files is provided below.

The model outputs are organised in a consistent way for each model, under
`./model_outputs/MODEL_NAME` (currently just `./model_outputs/kalbarri_2_coralbay`). 
Raster products typically include many tif files, but for each we also provide
a vrt file that can be used in GDAL based open-source GIS (e.g. QGIS, R, ...)
and will treat the tiles in unison.
* `arrival_time` contains arrival time rasters for each source zone, including the minimum arrival time and the scenario average arrival time.
  * For any scenario, the arrival time is defined as the time that the model first exceeds 1 cm above the background sea level. At any particular site, some scenarios may never meet this criteria. The minimum is derived from the minimum of all scenarios that do record an arrival time, and similarly for the average.
* `domains_shapefile` contains a shapefile with boxes showing the model domains and their resolutions.
* `elevation_in_model_no_tidal_adjustment` contains the elevation in AHD
* `elevation_in_model_with_tidal_adjustment` contains the elevation, adjusted so that 0m is a local high tide.
* `elevation_source_file_index` contains information on the datasets used to make the elevation in different areas.
* `tidal_adjustment` contains the tidal adjustment used in the modelling, equivalent to assuming the tsunami occurs at a local high tide.
* `inundation_rate_logic_tree_mean` contains the modelled rate of inundation (events / year) at the logic-tree-mean.
* `inundation_rate_84pc` contains the modelled rate of inundation (events / year) at the 84th percentile epistemic uncertainty. For computational efficiency this is not provided on the global domains.
* `inundation_rate_16pc` contains the modelled rate of inundation (events / year) at the 16th percentile epistemic uncertainty. For computational efficiency this is not provided on the global domains.
* `jatwc_inundation_zones` contains the modelled no-threat, marine warning and land-warning zones for ATWS zones resolved by the model, and the marine-warning maximum stage (in AHD).
* `max_depth_1in2500_84pc` contains the maximum depth (m) with exceedance-rate of 1/2500 (events/year) at the 84th percentile epistemic uncertainty. For computational efficiency results are only provided relatively near to the area of interest.
* `max_stage_AHD_1in2500_84pc` contains the maximum stage (m AHD) with exceedance-rate of 1/2500 (events/year) at the 84th percentile epistemic uncertainty. For computational efficiency results are only provided relatively near to the area of interest.
* `max_speed_1in2500_84pc` contains the maximum speed (m/s) with exceedance-rate of 1/2500 (events/year) at the 84th percentile epistemic uncertainty. For computational efficiency results are only provided relatively near to the area of interest.
* `max_flux_1in2500_84pc` contains the maximum flux (m^2/s) with exceedance-rate of 1/2500 (events/year) at the 84th percentile epistemic uncertainty. For computational efficiency results are only provided relatively near to the area of interest.

Separately we provide a folder `./ATWS_ZONES` which shows the tsunami warning zones used by JATWC.

## For more information

Code used for this project is in Geoscience Australia's ptha repository: `https://github.com/GeoscienceAustralia/ptha/tree/master/misc/WA_2025_2027`
* For Kalbarri to Coral Bay, see the subfolder `kalbarri_2_coralbay`
