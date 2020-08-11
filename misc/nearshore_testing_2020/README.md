This folder contains code and links to data required for a study comparing models and observations of nearshore tsunamis in Australia. The working title of the paper is *Davies, G. et al. (2020) Global dissipation models for simulating tsunamis at far-field coasts up to 60 hours post-earthquake: Multi-site tests in Australia.*

The sub-folders here are:

* [./analysis](./analysis) Analysis of the model results and creation of plots (this is the last step in the analysis) 
* [./breakwalls](./breakwalls) Make the breakwalls that we burn into our elevation model.
* [./elevation](./elevation) Elevation grids created for our model and used by the code in [./swals](./swals). The content of this folder should be augmented with data which can be downloaded as explained therein.
* [./sources](./sources) The source models used to initialise our tsunami model. The content of this folder should be augmented with data which can be downloaded as explained therein.
* [./swals](./swals) The tsunami model application specific code, job submission scripts, and some plotting scripts. 

The sub-folders contain code, and instructions to download zipped-versions of the folders with datasets (model input / model output / test data). The latter are too large to version control. In principle you could use this data and code to re-run the study. Beware our job-submission and module-loading scripts are specific to the computer we were using, and some of our R scripts included hard-coded links to our SWALS install on that machine - these kinds of things would have to be updated to run on another machine.

We rely on a separate install of the [SWALS source code](https://github.com/GeoscienceAustralia/ptha/tree/master/propagation/SWALS), and the [rptha](https://github.com/GeoscienceAustralia/ptha/tree/master/R) R package. Those links are also in the [ptha repository](https://github.com/GeoscienceAustralia/ptha). Various R packages from CRAN will also be used (generally available by calling `install.packages(...)` from within R).
