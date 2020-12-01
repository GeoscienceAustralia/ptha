# Tonga probabilistic tsunami inundation hazard
-----------------------------------------------

This folder contains code and data for an inundation PTHA for Tonga that was
initially developed in November 2020.

See README documentation within these folders for more information.

* [./analysis](./analysis) - Complex post-procressing of many-scenario hydrodynamic model outputs (e.g. probabilistic inundation calculations)
* [./elevation](./elevation) - Elevation data needed for hydrodynamic model. It includes various post-processing scripts to combined datasets into a format suitable for out hydrodynamic model.
* [./gauges](./gauges) - Gauge observations at Nuku'alofa, including de-tiding. Also locations of gauge-output points for hydrodynamic model.
* [./sources](./sources) - Earthquake-tsunami initial conditions for the model. This includes both: A) sources similar to historic events [for model testing], and; B) code to generate random tsunami scenarios for the probabilistic hazard assessment.
* [./swals](./swals) - Hydrodynamic model code, and basic post-processing code (e.g. creation of max-depth rasters; code to plot model-vs-observations at gauges; etc.)

