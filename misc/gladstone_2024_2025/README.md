# Gladstone 2024/25 Tsunami Inundation

The SWALS model in this project was developed for a probabilistic tsunami hazard assessment for the Gladstone region. The project was undertaken by Geoscience Australia in collaboration with Queensland Fire Department in 2024-25. The model places high-resolution grids in towns near Gladstone, from Baffle Creek to Yeppoon, and includes Heron Island, Lady Elliot Island and One Tree Island. An outline of how to setup the model is provided in [doc/general_guidance_on_model_setup.md](doc/general_guideance_on_model_setup.md). 

The code for this project was in large part modified from previous projects in [NSW](../nsw_2023_2024) and [WA](../SW_WA_2021_2024).
The scenario sampling method is more similar to the NSW study as it uses importance sampling based on offshore wave-heights in PTHA18.
However, it only has results from a single batch so doesn't require multiple importance sampling like in NSW (because in this case the model resolves inundation over a smaller coastline).

The key innovations in this study were a) the treatment of spatially varying highest astronomical [tide](tides) (HAT), which adjusts the elevation to account for the unevenness in HAT while ensuring the initial water level is quiescent. This is the subject of a conference paper in [Coasts and Ports 2025](https://coastsandports2025.com.au/). b) Nesting domains with a tree structure, rather than linear cascade. Customising the `nesting_level_parent_index` in [multidomain design mod](swals/model_multidomain_design_mod.f90) permits domains to nest inside any compatible parent domain. c) Using variable [friction](friction) based on land cover. 
The project uses [makefiles](makefile) to compile the input files for SWALS (e.g. breakwalls, elevation, friction, ...). It also includes single scenario analysis for the effect of 0.8 m of sea level rise and the presence/absence of friction from mangroves.


## Key folders

**Inputs to SWALS**
* [elevation](elevation) - Elevation data and preference order for the SWALS model
* [friction](friction) - Data and code used to define the spatially varying Manning roughness coefficient
* [multidomain_design](multidomain_design) - Create the boxes defining nested domains in the Gladstone region model's multidomain
* [breakwalls](breakwalls) - Create lines with elevation maxima that are burned into the elevation model, to enforce breakwalls and other local high-points irrespective of model resolution.
* [inverts](inverts) - Create lines with elevation minima that are burned into the model, to enforce inverts like drains
* [gauges](gauges) - Locations of point-gauge output used in the SWALS model
* [sources](sources) - Tsunami source models (water surface perturbations) for hazard scenarios, historical tsunamis, and other test cases.
* [initial_stage](initial_stage/) - Used to modify the model's initial water level to keep inland areas that are below sea level dry, if judged realistic to do so. 

**Running SWALS**
* [swals](swals) - Hydrodynamic model and outputs. Follow along in the [readme file](swals/README.md) to see the model development and workflow.

**Analysis of SWALS simulations**
* [analysis](analysis) - Computations that use multiple simulated hazard scenarios to produce outputs. The outputs include probabilistic hazard products, and inundation maps that correspond to JATWC categories.

## Modules
The project uses two sets of environment modules on NCI Gadi: [modules_R_431.sh](modules_R_431.sh) for working in R and with gdal; as well as [modules_SWALS_ifx_2024.sh](modules_SWALS_ifx_2024.sh) for compiling and running SWALS in fortran. Load them all as e.g.
```bash
source modules_R_431.sh
```
before running R scripts and likewise for fortran. The scripts may work for other versions of these software, but it's not guaranteed.

## Make SWALS inputs

The [makefile](makefile) in this directory is designed to make almost all the inputs for the SWALS (Shallow WAter Like Solvers) model. **Beware it does not account for everything (discussed below) so while helpful it is not a substitute for understanding what you've done**. Run `make` to check that the inputs that the script accounts for are up to date using GNU make. This requires the [SWALS modules](modules_SWALS_ifx_2024.sh) being loaded. Also try `make --dry-run` or `make --debug` first if you want to see what's not up to date before running. Using the makefile ensures that many of the hundreds of inputs are up to date, if they or their dependencies are modified. Inspecting the log file from `make -p > make.log` reveals there's about 73 csv files, 5 txt file, 54 shp files, and 35 R files that `make` tracks.

Beware that the makefile does not check every input. For example, many rasters are generated in GIS software. These rasters are [listed](elevation/make_swals_elevation_files_preference_list.R) in another [text file](elevation/swals_elevation_files_in_preference_order.txt) for SWALS. `Make` only sees changes in the text file, not the rasters themselves. Also, it currently doesn't include [source generation](sources/hazard/create_initial_conditions_for_scenarios.R) because this may require non-trivial compute resources. It may also require compute to [generate the friction rasters](friction/make_friction_rasters.R) (if required to regenerate).
