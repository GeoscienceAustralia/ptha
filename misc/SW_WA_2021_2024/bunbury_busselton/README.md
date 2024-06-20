Tsunami modelling for Busselton and Bunbury.

Updated variants of these codes were developed [in subsequent revisions](../greater_perth_revised2023). The latter are improved and simplified, so likely provide a better starting point for further adaptation.

Key folders are:
* [./analysis](./analysis) - Computations that use multiple simulated scenarios to produce outputs. Such as hazard products, and inundation maps that correspond to JATWC categories. **This was not ultimately used due to model updates**
* [./analysis_shutfloodgate](./analysis_shutfloodgate) - Computations that use multiple simulated scenarios to produce outputs. Such as hazard products, and inundation maps that correspond to JATWC categories. **This was used to represent outputs at Bunbury with a shut floodgate, but not at Busselton because there we later got updated data**
* [./analysis_NewVasseDrainOpenBunburyFloodgate](./analysis_NewVasseDrainOpenBunburyFloodgate) - Computations that use multiple simulated scenarios to produce outputs. Such as hazard products, and inundation maps that correspond to JATWC categories. **This was used to represent outputs at Bunbury with an open floodgate, and at Busselton (with the most up to date data there)**
* [./breakwalls](./breakwalls) - Create 3d lines that are burned into the elevation model, to enforce breakwalls and other local high-points irrespective of model resolution. 
* [./elevation](./elevation) - Elevation data + preference order for the SWALS model
* [./gauges](./gauges) - Locations of point-gauge output used in the SWALS model
* [./multidomain_design](./multidomain_design) - Create the boxes defining domains in the multidomain
* [./sources](./sources) - Tsunami source models
* [./swals](./swals) - Hydrodynamic model and outputs
