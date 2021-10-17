Code for a paper with working title: "From offshore to onshore PTHA: Efficient Monte-Carlo sampling"
-----------------------------------------------------------------------------------------------------

The directories contain:

* [./elevation](./gauges) -- Elevation data used for modelling.

* [./gauges](./gauges) -- Tide-gauge data used for model testing, and gauge locations stored by model.

* [./optimal_sampling](./optimal_sampling) -- A script that implements the offshore analysis of various Monte-Carlo schemes (prior to high-resolution inundation simulation). It creates Figures 2,3,4,5,6,7,8 from the paper (and many others that are not used). It also reports various statistics.

* [./sources](./sources) -- Code to select the particular set of random scenarios that was used for the onshore PTHA, and information on the source-models used for comparison with historic events. 

* [./swals](./swals) -- The tsunami model application code. This can be run after all the random source models have been created [using the code here](./sources/random/). As well as making simulations for all random scenarios, the code in this folder produces plots of the model domain and the comparison with historical events, which were manually mosaiced to make Figure 9 of the paper.

* [./analysis](./analysis) -- Probabilistic computations using the high-resolution model outputs. This can be run after all the simulations in the [./swals](./swals) directory are completed, and make Figures 10, 11, 12 in the paper.

* [./regional_setting](./regional_setting) -- Code to make panels of Figure 1 in the paper.
