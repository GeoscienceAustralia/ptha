# Gladstone 2024/25 Tsunami Inundation

This model is being developed for to develop probabilistic inundation hazards for the Gladstone region for Queensland Fire Department in 2024-25. This model places high-resolution grids in towns near Gladstone, from Baffle Creek to Yeppoon and includes several offshore islands. An outline of how to setup the model is provided in [./doc/general_guidance_on_model_setup.md](./doc/general_guidance_on_model_setup.md).

## Key folders

**Inputs to SWALS**
* [elevation](./elevation) - Elevation data + preference order for the SWALS model
* [friction](/friction) - Data containing the extent of mangroves
* [multidomain_design](./multidomain_design) - Create the boxes defining domains in the multidomain
* [breakwalls](./breakwalls) - Create lines that are burned into the elevation model, to enforce breakwalls and other local high-points irrespective of model resolution.
* [inverts](./inverts) - Lines with elevation burned into the model to enforce inverts like drains
* [gauges](./gauges) - Locations of point-gauge output used in the SWALS model
* [sources](./sources) - Tsunami source models
* [initial_stage](./initial_stage/) - Modified initial water level for dry inland areas below sea level 

**Running SWALS**
* [swals](./swals) - Hydrodynamic model and outputs. Follow along in the [readme file](swals/README.md) to see the model development and workflow.

**Analysis of SWALS simulations**
* [analysis](./analysis) - Computations that use multiple simulated scenarios to produce outputs. Such as probabilistic hazard products, and inundation maps that correspond to JATWC categories.

## Make SWALS inputs

The [makefile](./makefile) in this directory is designed to make all the inputs for the SWALS (Shallow WAer Like) model. Run `make` to ensure all inputs are up to date using GNU make. Also try `make --dry-run` or `make --debug` first if you want to see what's not up to date before running. It currently doesn't include source generation because this may require non-trivial compute resources. It may also require compute to generate the friction rasters (if required to regenerate).

Using the makefile ensures that many of the hundreds of inputs are up to date, if they, or their dependencies are modified. Inspecting the log file from `make -p > make.log` reveals there's about 73 csv files, 5 txt file, 54 shp files, and 35 R files that make tracks. However, beware that it does not compile every input. For example, many rasters are generated in GIS software, and then another file might just point to its file location. `Make` only sees changes in the file with the list.
