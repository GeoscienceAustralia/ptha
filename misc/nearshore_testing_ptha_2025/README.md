# Code and data for a study that tests PTHA18 stochastic earthquake-tsunami models at many tide gauges in Australia.
--------------------------------------------------------------------------------------------------------------------

This folder contains code and links to data used for the paper INSERT-NAME-AND-LINK-WHEN-PUBLISHED.

The information is provided for transparency and to assist future studies.
It was partly run on the NCI Gadi supercomputer, and partly on the authors
local machine. The code makes assumptions about those environments that would
need to be changed if you were trying to run the code elsewhere (mainly the
location of installed software, and the location of datasets that existed at a
higher level in the filesystem during the analysis, as compared to where it is
provided here).

## Subfolders
* `./analysis_ptha18_scenarios_2025` - Analysis of the tsunami models, which can be run after all calculations in the `swals` folder are complete.
* `./breakwalls` - Information on small scale linear features that are burned into the model elevation, such as breakwalls that would not otherwise be captured at the SWALS model's native resolution.
* `./elevation` - Links to input elevation data used by the model, and also the elevation as seen by the SWALS model (derived from the input elevation data using the model code in `./swals`)
* `./gauges` - Tide gauge data and a script that provides a uniform interface. Beware that scripts in this project use a different file path to refer to the tide gauge interface script (since when the calculations were implemented, this tide gauge data was stored in a central location that cannot easily be represented in this repository). 
* `./ptha18_scenarios_random` - Link to code used to sample random scenarios from PTHA18, as well as the actual scenarios.
* `./sources` - Sources for historical tsunamis derived from published source inversions
* `./swals` - Code used to model the tsunami for all scenarios
