# Combine inundation rates (84th percentile model) from multiple models into a single raster output.

This code creates a single inundation-rate (84th percentile) raster using the "best" model at each site.
In some areas we use the maxima of multiple models.
- Bunbury behind the storm surge barrier: Maximum of models with floodgate open/shut
- Geraldton: Maximum of Midwest model (including buildings in topography) and Greater Perth model (excluding buildings in topography).

It also creates a contour with the 1/2500, 1/500 and 1/100 inundation rate (84th percentile). 

## Folders

* merged_inundation_rate_84pc -- contains output rasters with the merged inundation rate (84th percentile), including a vrt that can can be used to view them all at once.
* merged_inundation_rate_contour_1in2500_1in500_1in100_84pc -- A shapefile that contours the aforementioned rasters. The contours were created in parallel and so are disconnected at raster boundaries.
* merged_model_raster_output_zones -- Polygons defining the locations and cell-size of output rasters
* priority_model_zones -- User-defined polygons used to tell the code which model should be used in which area.

## How to run

First create polygons defining regions where rasters will be created, in part using shapefiles in `priority_model_zones`.

    Rscript define_raster_output_regions.R

Then make the rasters and contour (this uses parallel computing with hard-coded mc.cores)

    Rscript make_rasters.R
