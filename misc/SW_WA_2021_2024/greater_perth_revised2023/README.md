Tsunami modelling for the Greater Perth Area, 2023 update (part of the WA
Tsunami Inundation Modelling Project with DFES).

This folder was used to re-run some previous [Greater Perth
models](../greater_perth_metro_area) after we obtained improved elevation
data for a few sites of interest. The results are almost identical to the previous
models at most sites, but have small changes near the updated elevation. 

But the code here has been substantially improved and simplified to make additional adaptation easier.

Additional useful content may be found in the [Greater Perth Modelling
files](../greater_perth_metro_area). For instance that includes checks of
sea-levels and initial model developments that are also useful here.

Key folders are:
* [./analysis](./analysis) - Computations that use multiple simulated scenarios to produce outputs. Such as hazard products, and inundation maps that correspond to JATWC categories.
* [./breakwalls](./breakwalls) - Create 3d lines that are burned into the elevation model, to enforce breakwalls and other local high-points irrespective of model resolution. 
* [./elevation](./elevation) - Elevation data + preference order for the SWALS model
* [./inverts](./inverts) - Create 3d lines that are burned into the elevation model, to enforce channels and other local low-points irrespective of model resolution. 
* [./gauges](./gauges) - Locations of point-gauge output used in the SWALS model
* [./multidomain_design](./multidomain_design) - Create the boxes defining domains in the multidomain
* [./sources](./sources) - Tsunami source models
* [./swals](./swals) - Hydrodynamic model and outputs

