# Definitive outputs for the WA tsunami inundation modelling project, 2021-2024

This folder contains code to collate outputs from multiple tsunami models
developed during the project. The intention is to make it easier for future
users to find the most up-to-date outputs of interest.

The modelling was conducted in three phases, each with particular focus areas:
* Phase 1, Greater Perth, see `./model_outputs/greater_perth`
  * Two Rocks to Dawesville including Rottnest Island and Garden Island
  * Bunker Bay (near Cape Naturaliste)
* Phase 2, Bunbury/Busselton, see `./model_outputs/bunbury_busselton`
  * Bunbury, Busselton and Dunsborough
  * This model has the Bunbury floodgate open. The Busselton floodgates are always closed.
    * A separate model provides results with Bunbury floodgate closed. This variant should not be used outside of Bunbury, and in parts of Busselton it uses out-of-date elevation data. 
      * See `./model_outputs/bunbury_busselton_alternative_with_bunbury_floodgate_shut`
* Phase 3, Midwest, see `./model_outputs/midwest`
  * Towns from Guilderton to Geraldton, and the Abrolhos Islands.
  * At Geraldton, the Midwest model includes a few buildings in the elevation data (to investigate how they may affect inundation). 
    * If you don't want this, then Geraldton results without buildings are available in the Phase 1 and Phase 2 models (high resolution). 

The three models differ in the areas treated at high resolution, and the
underlying datasets (to a degree). There is also lots of overlap between
the three models, and their results are typically quite consistent.

We have also combined the inundation rates (84th percentile) from multiple
models into a single output product. 
* See `./combined_models`, with subfolders including 
  * The inundation rates (84th percentile) in raster format in `./combined_models/merged_inundation_rate_84pc`
  * Contours of the inundation rates (84th percentile) at 1/100, 1/500 and 1/2500 in line shapefile format in `./combined_models/merged_inundation_rate_contour_1in2500_1in500_1in100_84pc`
* For details on which model is used where, see `./combined_models/README.md`.

## Structure of the model output folders

All files are loaded in a QGIS session contained in the folder `./QGIS_session`. The same folder contains a backup copy of the QGIS session.

If the QGIS session doesn't work for you, or you'd rather use them in another way, information on the files is provided below.

The model outputs are organised in a consistent way for each model, under `./model_outputs/MODEL_NAME` (e.g. for Greater Perth it is `./model_outputs/greater_perth`). Raster products typically include many tif files, but for each we also provide a vrt file that can be used in GDAL based open-source GIS (e.g. QGIS, R, ...) and will treat the tiles in unison.
* `arrival_time` contains arrival time rasters for each source zone, including the minimum arrival time and the scenario average arrival time.
  * For any scenario, the arrival time is defined as the time that the model first exceeds 0.61 m (i.e. 1 cm above the background sea level for almost all of the model domain). At any particular site, some scenarios will never meet this criteria. The minimum is derived from the minimum of all scenarios that do record an arrival time, and similarly for the average.
* `domains_shapefile` contains a shapefile with boxes showing the model domains and their resolutions.
* `elevation_in_model` contains the elevation as seen by the model
* `inundation_rate_logic_tree_mean` contains the modelled rate of inundation (events / year) at the logic-tree-mean.
* `inundation_rate_84pc` contains the modelled rate of inundation (events / year) at the 84th percentile epistemic uncertainty. For computational efficiency this is only provided for relatively high resolution domains.
* `inundation_rate_16pc` contains the modelled rate of inundation (events / year) at the 16th percentile epistemic uncertainty. For computational efficiency this is only provided for relatively high resolution domains.
* `jatwc_inundation_zones` contains the modelled no-threat, marine warning and land-warning zones for ATWS zones resolved by the model.
* `max_depth_1in2500_84pc` contains the maximum depth (m) with exceedance-rate of 1/2500 (events/year) at the 84th percentile epistemic uncertainty. For computational efficiency results are only provided for tiles with good resolution, and the min/max are truncated.
* `max_stage_1in2500_84pc` contains the maximum stage (m) with exceedance-rate of 1/2500 (events/year) at the 84th percentile epistemic uncertainty. For computational efficiency results are only provided for tiles with good resolution, and the min/max are truncated.
* `max_speed_1in2500_84pc` contains the maximum speed (m/s) with exceedance-rate of 1/2500 (events/year) at the 84th percentile epistemic uncertainty. For computational efficiency results are only provided for tiles with good resolution, and the min/max are truncated.
* `max_flux_1in2500_84pc` contains the maximum flux (m^2/s) with exceedance-rate of 1/2500 (events/year) at the 84th percentile epistemic uncertainty. For computational efficiency results are only provided for tiles with good resolution, and the min/max are truncated.

Separately we provide a folder `./ATWS_ZONES` which shows the tsunami warning zones used by JATWC.

## For more information

See the report to understand the products (copy provided in this folder).

Code used for this project is in Geoscience Australia's ptha repository: 

`https://github.com/GeoscienceAustralia/ptha/tree/master/misc/SW_WA_2021_2024`
  * For Greater Perth, see the subfolder `greater_perth_revised2023`
  * For Bunbury Busselton, see the subfolder `bunbury_busselton`
  * For Midwest, see the subfolder `guilderton_to_geraldton`
