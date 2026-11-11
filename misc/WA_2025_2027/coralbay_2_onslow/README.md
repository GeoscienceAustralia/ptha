# Inundation PTHA for Coral Bay to Onslow

Part 2 of phase 4 of the Western Australia Tsunami Inundation Modelling Project with DFES.

The project sets up the inputs (inverts, breakwalls, gauges, elevation and multidomains) for SWALS.
Once the inputs are ready, compile SWALS with it's [makefile](./swals/makefile) and run from there. 
The project uses [one module file](./modules_R_431.sh) to source for running R scripts,
and another [one](./modules_ifx_2025_llvm.sh) to compile SWALS.

Most of the code was derived from the first half of the modelling [../kalbarri_2_coralbay](../kalbarri_2_coralbay/) to cover the new region. The code was reorganised in several places:
- test and validation PBS scripts in [swals](./swals) were placed into the test and historic subdirectories. 
- The gauge plots comparing model recreations to historic records have been moved from within [swals](./swals) to [analysis/gauge_plots](./analysis/gauge_plots).
- The project has only two module files in the root directory: [./modules_ifx_2025_llvm.sh](./modules_ifx_2025_llvm.sh) for compiling swals and [./modules_R_431.sh](./modules_R_431.sh) for postprocessing. The R modules file now include `imagemagick` and `ffmpeg` for making animations.
- The source directory is a symlink to the same sources used in [../kalbarri_2_coralbay](../kalbarri_2_coralbay/).
- The namelists in [./swals/multidomain_design_control](./swals/multidomain_design_control) have simplified names.
- PBS scripts have the file extension `.pbs` (as opposed to `.sh`) and their filenames may have shortened to remove a leading `run_` or `compute_`.   

New features for this project are:
- the ability to create an [animation](./animation/),
- the ability to export elevation layers from QGIS into a file directly for SWALS using [./elevation/save_qgis_layers.py](./elevation/save_qgis_layers.py) from within QGIS.
- explicit testing scripts to initialise the load balance process in [swals/test](./swals/test). Start with `small_TEST_ONLY_initialise.pbs` then `small_TEST_ONLY_check_load_balance.pbs` to get a load balance, possibly with iteration.
- an additional test [./swals/test/instability_carnarvon.pbs](./swals/test/instability_carnarvon.pbs) based on an instability in the first half of the modelling.
