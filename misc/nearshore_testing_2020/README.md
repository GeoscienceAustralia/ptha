Testing various alternative tsunami modelling approaches at Australian tidal gauges, using published source-inversions to represent the historic tsunami sources.
--------------------------------------------------------------------------------------------------------

This folder contains code and links to data required for a study comparing models and observations of nearshore tsunamis in Australia. The work was published in this paper: [Davies, G. et al. (2020) Global dissipation models for simulating tsunamis at far-field coasts up to 60 hours post-earthquake: Multi-site tests in Australia.](https://www.frontiersin.org/articles/10.3389/feart.2020.598235/full).

The sub-folders contain code, and instructions to download zipped-versions of the folders with datasets (model input / model output / test data). The latter are too large to version control. The folders are:

* [./breakwalls](./breakwalls) Make the breakwalls that we burn into our elevation model.
* [./gauges](./gauges) Processed tide-gauge data used in the study, along with a script which packs it all into a convenient data-structure (which is used throughout R scripts in the other folders). The content of this folder should be augmented with data which can be downloaded, as explained therein.
* [./elevation](./elevation) Elevation grids created for our model and used by the code in [./swals](./swals). The content of this folder should be augmented with data which can be downloaded, as explained therein.
* [./sources](./sources) The source models used to initialise our tsunami model. The content of this folder should be augmented with data which can be downloaded, as explained therein.
* [./swals](./swals) The tsunami model application specific code, job submission scripts, and some plotting scripts. 
* [./analysis_nearshore2020](./analysis_nearshore2020) Analysis of the model results and creation of most plots (this is the last step in the analysis). The content of this folder should be augmented with data which can be downloaded, as explained therein. This includes [figures comparing all models and gauges for all sources](http://dapds00.nci.org.au/thredds/fileServer/fj6/PTHA/Nearshore_testing_2020/all_models_sites_vs_data.zip). The README in this folder also includes pointers to code that was used to make each figure in the paper.


In principle you could use this data and code to re-run the study. In that instance, beware our job-submission and module-loading scripts are specific to [NCI's Gadi computer](https://nci.org.au/our-systems/hpc-systems), and might need to be changed for other machines using other queueing systems. Also some of our R scripts included hard-coded links to a copy of the [ptha](https://github.com/GeoscienceAustralia/ptha) repository on that machine - these kinds of things would have to be updated to run on another machine.

We rely on a separate install of the [rptha](https://github.com/GeoscienceAustralia/ptha/tree/master/R) R package, and repeatedly link to a copy of the [ptha repository](https://github.com/GeoscienceAustralia/ptha) to use the SWALS source code. Various R packages from CRAN are also used (and these can generally be installed by calling `install.packages(...)` from within R).

The code and data herein is provided for scientific transparency (i.e. to show what we did). We do not guarentee that the contents are totally accurate or complete, or that they are suitable for any purposes other than understanding what was done in our study. You should not solely rely on this information when making a commercial decision .
