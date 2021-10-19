Code for a paper with working title: "From offshore to onshore PTHA: Efficient Monte-Carlo sampling"
-----------------------------------------------------------------------------------------------------

This folder contains application code used for our manuscript on Monte-Carlo sampling for PTHA. The directories are:

* [./elevation](./gauges) -- Intended to contain elevation data used for modelling (rasters not included here).

* [./gauges](./gauges) -- Tide-gauge data used for model testing, and gauge locations stored by the model.

* [./optimal_sampling](./optimal_sampling) -- A script that implements the offshore analysis of various Monte-Carlo schemes (prior to high-resolution inundation simulation). It creates Figures 2,3,4,5,6,7,8 from the paper (and many others that are not used in the paper). It also reports various statistics.

* [./sources](./sources) -- Code to select the particular set of random scenarios that were used for the onshore PTHA computations. Also information on the source-models used for comparison with historic events. 

* [./swals](./swals) -- The tsunami model application code. This can be run after all the random source models have been created [using the code here](./sources/random/). As well as making simulations for all random scenarios, the code in this folder produces plots of the model domain and the comparison with historical events, which were manually mosaiced to make Figure 9 of the paper.

* [./analysis](./analysis) -- Probabilistic computations using the high-resolution model outputs. This can be run after all the simulations in the [./swals](./swals) directory are completed, and make Figures 10, 11, 12 in the paper.

* [./regional_setting](./regional_setting) -- Code to make panels of Figure 1 in the paper.


While intending to make the analysis transparent, the code is application-specific and non-trivial edits would be required to run on other machines, or adapt to other problems. The code was mostly run on [NCI's Gadi computer](https://nci.org.au/our-systems/hpc-systems), and our job-submission and module-loading scripts are specific to that machine. Some of our R scripts also include hard-coded links to data that exists on the machines we used, but cannot be provided here. 
