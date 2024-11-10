# NSW tsunami modelling project results

This folder contains outputs related to the NSW tsunami inundation modelling
project.

## Contents

* `NSW_Notes_on_modelling_2024_11_06.pdf` contains PPT slides from a project meeting. 
  * PURPOSE: Help understand the background and motivation for the different products provided here.

* `QGIS_session/` includes a QGIS project that loads all of the datasets below.
  * You don't have to use this but it might be convenient.

* `domains_shapefile/` is a shapefile showing the hydrodynamic model "domains". 
  * PURPOSE: Quick overview of the modelled regions and their resolutions.
  * Each domain is a rectangular region in which we model the tsunami with a given spatial resolution. 
    * The model includes hundreds of domains which communicate with each other as the tsunami evolves (two-way nesting). 
    * In combination these domains enable the model to have near-global coverage, and to resolve inundation (~ 35m cells) for the entire NSW coast and Lord Howe Island.
  * Our raster outputs typically include hundreds of tifs (one per domain for each variable). 
    * `domains_shapefile` can help you quickly understand where different tifs are located, if necessary.

* `elevation_in_model/` contains elevation rasters derived from the model.
  * PURPOSE: Examine the the elevation according to the model.
  * There is one tif per domain
    * The tif boundaries are slightly outside the rectangular boundaries in `domains_shapefile` (reflecting halos used for parallel communication in our model) but there is no ambiguity because in the regions of overlap, only 1 tif will contain non-missing data.
  * There is also a single VRT file, which seamlessly treats the set of tifs like a single file. 
    * It was created via GDAL with a command similar to `gdalbuilvrt -resolution highest NAME_OF_OUTPUT_VRT_FILE.vrt elevation0*.tif`
    * This can be used in most open-source GIS software, and works well in my experience with QGIS, R, and GDAL. 
    * For ESRI software you might have more success by creating a "Raster Mosaic dataset" with the tifs (albeit I haven't tried this).

* `JATWC_inundation_zones/` contains modelling products that are linked to JATWC tsunami warning category.
  * PURPOSE: Modelling products that relate our inundation results to JATWC warnings. They will be a good starting point for the design of evacuation zones.
  * For further information see the README inside that folder, and the presentation in the current directory.

* `Probabilistic_inundation_hazard/` contains modelling products that depict:
  * The rate of inundation (events/year)
  * Flow variables at a specified rate of inundation 
    * e.g. `depth above initial condition` with a 1/2500 chance of exceedance at the 84th percentile epistemic uncertainty
  * Flood hazard categories
  * PURPOSE: These products are a key input used to create the products inside `JATWC_inundation_zones/`. They also provide context that could help with the design of evacuation zones.
  * For further information see the README inside that folder, and the presentation in the current directory.

* `Arrival_times/` contains modelled minimum and average arrival times:
  * The arrival times are in seconds post-earthquake
  * Tsunami arrival is defined as the first time at which the modelled tsunami exceeds 1cm above the background sea level of 1.1m AHD. Note this can lead to discontinuous arrival times.
  * The mean/minimum arrival time is the average/minimum respectively over all modelled scenarios from the source zone, ignoring their likelihoods.
