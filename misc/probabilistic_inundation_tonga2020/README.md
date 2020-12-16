# Tonga probabilistic tsunami inundation hazard
-----------------------------------------------

This folder contains code and data for an inundation PTHA for Tonga that was
initially developed in November 2020.

Rasters depicting the tsunami inundation depth around Tongatpu with a 10% and 2% chance of exceedance in 50 years are provided for download. We provide the [results using a background sea-level of 0m MSL](http://dapds00.nci.org.au/thredds/fileServer/fj6/PTHA/Tonga_2020/alternate_ptha18_tonga_MSL0.zip), as well as the same [results using a more conservative background sea-level of 0.8m MSL](http://dapds00.nci.org.au/thredds/fileServer/fj6/PTHA/Tonga_2020/alternate_ptha18_tonga_MSL0.8.zip). Some additional information is provided in the [probabilistic_inundation](./analysis/probabilistic_inundation) code and documentation.

# Background on these files

The codes were setup to be run on NCI's Gadi machine, and while they can certainly
be adapted to run elsewhere, some parts of the code include hard-coded links to
other files (e.g. the location of the ptha repository on Gadi) that would need to be updated.

See README documentation within these folders for more information.

* [./analysis](./analysis) - Complex post-procressing of many-scenario hydrodynamic model outputs (e.g. probabilistic inundation calculations)
* [./elevation](./elevation) - Elevation data needed for hydrodynamic model. It includes various post-processing scripts to combined datasets into a format suitable for out hydrodynamic model.
* [./gauges](./gauges) - Gauge observations at Nuku'alofa, including de-tiding. Also locations of gauge-output points for hydrodynamic model.
* [./sources](./sources) - Earthquake-tsunami initial conditions for the model. This includes both: A) sources similar to historic events [for model testing], and; B) code to generate random tsunami scenarios for the probabilistic hazard assessment.
* [./swals](./swals) - Hydrodynamic model code, and basic post-processing code (e.g. creation of max-depth rasters; code to plot model-vs-observations at gauges; etc.)

As part of this study we compared the model with historical observations at Nuku'alofa tide-gauge. Here is one example, for a small tsunami following a local earthquake in 2006; the scenario is from PTHA18 but was previously found to show broad agreement with DART buoy results for this tsunami. Other examples, for tsunamis originating in Japan and Chile, are available [at the bottom of this page](./swals/).

![Model-vs-data plot for the Tonga 2006 earthquake-tsunami](swals/plots/historic_events_time_series_plots/Tonga2006/nukualofa_gauge_modelVdata_Tonga2006_validation_PTHA18_VAUS_26849_load_balance-risetime_0-ambientsealevel_0.0-full-linear_with_manning-0.035-highres_tonga-RUN_20201130_172740871.png)
