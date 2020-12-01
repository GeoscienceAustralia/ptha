# Elevation data -- raw and derived -- for use with the SWALS Tonga model
-------------------------------------------------------------------------

In this folder we collate and process elevation data to support the SWALS model for Tonga. The folder contains:

* [./for_model](./for_model) A set of rasters that will be used by the SWALS model to set the elevation (they are given a preference order in the code)
* [./walls](./walls) Code and data to compute elevation maxima along manually defined lines. These will be used in the SWALS model, where we enforce the elevation maxima along narrow features that may act as topographic barriers, irrespective of the grid size. 
