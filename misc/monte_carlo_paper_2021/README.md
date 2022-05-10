Code for the paper "[Davies et al., (2022) From offshore to onshore PTHA via efficient Monte Carlo sampling, Geophysical Joournal International](https://doi.org/10.1093/gji/ggac140)"
---------------------------------------------------------------------------------------

This folder contains application code used for our manuscript on Monte-Carlo sampling for PTHA. It is released to show how the calculations in the paper were implemented (i.e. for transparency), and is not expected to be developed or maintained in the longer term.

Note that a separate tutorial on using these sampling methods with PTHA18 is available [here](../../ptha_access/example_event_access_scripts/random_scenarios_non_uniform_and_importance_sampling).

The directories are:

* [./optimal_sampling](./optimal_sampling) -- A script that implements the offshore analysis of various Monte-Carlo schemes (prior to high-resolution inundation simulation). It creates Figures 2,3,4,5,6,7,8 from the paper (and many others that are not used in the paper). It also reports various statistics which summarise the performance of different Monte-Carlo schemes, and computes the non-uniform sampling effort within magnitude-bins that is used in the paper. 

* [./elevation](./elevation) -- Intended to contain elevation data used for hydrodynamic modelling (rasters not included here). To save space we provide a link to another part of the repository that has the same contents.

* [./gauges](./gauges) -- Tide-gauge data used for hydrodynamic model testing, and gauge locations stored by the model. 

* [./sources](./sources) -- Code to select the particular set of random scenarios that were used for the onshore PTHA computations. Also information on the source-models used for comparison with historic events.

* [./swals](./swals) -- The tsunami model application code, which uses the [SWALS solver](../../propagation/SWALS). The models here can be run after all the random source models have been created [using the code here](./sources/random/). As well as making simulations for all random scenarios, the code in this folder produces plots of the model domain and the comparison with historical events, which were manually combined to make Figure 9 of the paper.

* [./analysis](./analysis) -- Probabilistic computations using the high-resolution model outputs. This can be run after all the simulations in the [./swals](./swals) directory have completed. These codes make Figures 10, 11, 12 of the paper.

* [./regional_setting](./regional_setting) -- Code to make panels of Figure 1 in the paper (the panels were manually combined to produce Figure 1).

While intending to make the analysis transparent, the code here is application-specific. Non-trivial edits would be required to run on other machines, or adapt to other problems. The code was mostly run on [NCI's Gadi computer](https://nci.org.au/our-systems/hpc-systems), and our job-submission and module-loading scripts are specific to that machine. Some of our R scripts also include hard-coded links to data that exists on the machines we used, but cannot be provided here. 
