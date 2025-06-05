# Gladstone 2024/25 Tsunami Inundation

This model is being developed for to develop probabilistic inundation hazards for the Gladstone region for Queensland Fire Department in 2024-25. This model places high-resolution grids in towns near Gladstone, from Baffle Creek to Yeppoon and includes several offshore islands. An outline of how to setup the model is provided in [doc/general_guidance_on_model_setup.md](doc/general_guideance_on_model_setup.md).

The code for this project was in large part modified from previous projects in [NSW](../nsw_2023_2024) and [WA](../SW_WA_2021_2024).
The scenario sampling method is more similar to the NSW study as it uses importance sampling based on offshore wave-heights in PTHA18.
However, it only has results from a single batch so doesn't require multiple importance sampling like in NSW.

The key innovations in this study were a) the treatment of spatially varying highest astronomical [tide](tides) (HAT), which adjusts the elevation to account for the unevenness in HAT while ensuring the initial water level is quiescent. This is the subject of a conference paper in [Coasts and Ports 2025](https://coastsandports2025.com.au/). b) Nesting domains with a [tree structure](swals/model_multidomain_design_mod.f90), rather than linear cascade. c) Using variable [friction](friction) based on land cover. 
The project uses [makefiles](makefile) to compile the input files for SWALS (e.g. breakwalls, elevation, friction, ...). It also includes single scenario analysis for the effect of 0.8 m of sea level rise and the effect of the presence/absence of friction from mangroves. 


## Key folders

**Inputs to SWALS**
* [elevation](elevation) - Elevation data and preference order for the SWALS model
* [friction](friction) - Data containing the extent of mangroves
* [multidomain_design](multidomain_design) - Create the boxes defining domains in the multidomain
* [breakwalls](breakwalls) - Create lines that are burned into the elevation model, to enforce breakwalls and other local high-points irrespective of model resolution.
* [inverts](inverts) - Lines with elevation burned into the model to enforce inverts like drains
* [gauges](gauges) - Locations of point-gauge output used in the SWALS model
* [sources](sources) - Tsunami source models
* [initial_stage](initial_stage/) - Modified initial water level for dry inland areas below sea level 

**Running SWALS**
* [swals](swals) - Hydrodynamic model and outputs. Follow along in the [readme file](swals/README.md) to see the model development and workflow.

**Analysis of SWALS simulations**
* [analysis](analysis) - Computations that use multiple simulated scenarios to produce outputs. Such as probabilistic hazard products, and inundation maps that correspond to JATWC categories.

## Make SWALS inputs

The [makefile](makefile) in this directory is designed to make almost all the inputs for the SWALS (Shallow WAer Like) model. Run `make` to ensure all inputs are up to date using GNU make. Also try `make --dry-run` or `make --debug` first if you want to see what's not up to date before running. Using the makefile ensures that many of the hundreds of inputs are up to date, if they, or their dependencies are modified. Inspecting the log file from `make -p > make.log` reveals there's about 73 csv files, 5 txt file, 54 shp files, and 35 R files that `make` tracks.


Beware that it does not compile every input. For example, many rasters are generated in GIS software. These rasters are [listed](elevation/make_swals_elevation_files_preference_list.R) in another [text file](elevation/swals_elevation_files_in_preference_order.txt) for SWALS. `Make` only sees changes in the text file with the list. Also, it currently doesn't include [source generation](sources/hazard/create_initial_conditions_for_scenarios.R) because this may require non-trivial compute resources. It may also require compute to [generate the friction rasters](friction/make_friction_rasters.R) (if required to regenerate).
