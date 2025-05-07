# Code and data for a study that tests PTHA18 stochastic earthquake-tsunami models at many tide gauges in Australia.
--------------------------------------------------------------------------------------------------------------------

This folder contains code and links to data used for the paper INSERT-NAME-AND-LINK-WHEN-PUBLISHED.

The information is provided for transparency and to assist future studies.
However, the code was not designed to be naively re-run on other machines. It
was partly run on the NCI Gadi supercomputer, and partly on the authors
local machine. The codes make assumptions about those environments (e.g. the
location of data and installed software) that would need to be changed if you
were trying to run the code elsewhere. 

## Subfolders
* `./breakwalls` - Information on small scale linear features that are burned into the model elevation, such as breakwalls that would not otherwise be captured at the SWALS model's native resolution.
* `./elevation` - Links to input elevation data used by the model, and also the elevation as seen by the SWALS model (derived from the input elevation data using the model code in `./swals`)
* `./ptha18_scenarios_random` - Link to code used to sample random scenarios from PTHA18, as well as the actual scenarios.
* `./swals` - Code used to model the tsunami for all scenarios
