Tsunami modelling codes for Greater Perth Metro Area 
----------------------------------------------------

This folder contains code to run the tsunami models, including:
* Scenarios like historical events [(see here)](../sources/like_historic/)
* Random PTHA18 scenarios [(see here)](../sources/hazard/)

The key codes are described below; consult comments in the code itself for further documentation:

## Hydrodynamic model setup and compilation

* [model.f90](model.f90) and [model_local_routines.f90](model_local_routines.f90) are used to setup the hydrodynamic model.
    * [make_model_ifort](make_model_ifort) was used to compile the model on NCI, after loading the modules in [SWALS_ifort_modules_2022.sh](SWALS_ifort_modules_2022.sh) 
        * `source SWALS_ifort_modules_2022.sh; make -B -f make_model_ifort`
    * The model requires information on the multidomain layout ([created here](../multidomain_design/)) and the initial conditions. The latter include [the initial water-surface perturbation](../sources/) and [the elevation grids](../elevation/) and [the breakwall geometry](../breakwalls).

## Creation of qsub scripts to run the hydrodynamic model for random PTHA scenarios

* [create_random_ptha_qsub_scripts_sealevel60cm_reviseddomain_highres.R](create_random_ptha_qsub_scripts_sealevel60cm_reviseddomain_highres.R) is used to make qsub scripts which run the hydrodynamic model for all the random scenarios (for the full-resolution runs), and also `tar` the resulting output folders (to prevent creation of too many files on NCI). 
    * The script can be run with 
        * `Rscript create_random_ptha_qsub_scripts_sealevel60cm_reviseddomain_highres.R`.
    * It produces approximately 50 qsub scripts which can be separately submitted. One can also submit all the scripts at once using shell commands such as: 
        * `for i in run_greaterperth_sealevel60cm_reviseddomain_highres_*.sh; do qsub $i; mv $i revised_submitted_qsub_scripts; done`. 
* [create_random_ptha_qsub_scripts_sealevel60cm_lowres.R](create_random_ptha_qsub_scripts_sealevel60cm_lowres.R) is similar to above, but used to submit low-resolution versions of the models (for convergence testing). 
    * Before this is run, the model should be recompiled to use coarser resolution (using `mesh_refine=2_ip` in [model.f90](model.f90) ). 
* [DEBUG_run_with_old_nesting.sh](DEBUG_run_with_old_nesting.sh) is used to run a single random scenario which went unstable with the high-resolution run above. 
    * To do this the model was recompiled to use the old swals nesting scheme (compiling with `-DOLD_PROCESS_DATA_TO_SEND_B4FEB22`)
        * See commented out line in [make_model_ifort](make_model_ifort) for how to do this. 
    * From experience the different nesting schemes both work well in most cases and lead to unimportant changes in the results (this blow-up is a rare exception). 
    * The newer scheme has advantages when nesting fine staggered-grid models, which does not feature in this case.

## Creating rasters from the hydrodynamic model results

* [create_tarred_rasters_from_tarred_multidomains.R](create_tarred_rasters_from_tarred_multidomains.R) makes raster output files from the tarred multidomain directories, and puts them in a `tar` folder.
    * Note that GDAL can read rasters inside tar folders, [see here](https://gdal.org/user/virtual_file_systems.html)). 
    * The tarring is needed to avoid creating too many output files on NCI. 
    * Submit this script using [run_create_tarred_rasters_from_tarred_multidomains.R](run_create_tarred_rasters_from_tarred_multidomains.R), after editing the command to work on the appropriate folder.

## Simulating historic tsunamis

* [run_model_Fuji2004_6_nodes.sh](run_model_Fuji2004_6_nodes.sh) is used to submit a run similar to the 2004 Indian Ocean tsunami. 
    * See commented out variants of the commands therein for information on running a version with time-varying forcing (different for each unit source), a version with instantaneous forcing, and short versions used to create load balance files.
    * [plot_gauges_perth_sumatra2004.R](plot_gauges_perth_sumatra2004.R) is used to plot gauges from these runs. A post-processed gauge file has to have been created already using [plots/plot_sumatra2004.R](plots/plot_sumatra2004.R). 
* [run_Fuji_Sumatra2005_6_nodes.sh](run_Fuji_Sumatra2005_6_nodes.sh) is used to submit a run simular to the 2005 Sumatra tsunami. 
    * [plot_gauges_perth_sumatra2005.R](plot_gauges_perth_sumatra2005.R) is used to plot gauges from these runs. A post-processed gauge file has to have been created already using [plots/plot_sumatra2005.R](plots/plot_sumatra2005.R). 

## Other miscellaneous code

* [create_plots_from_tarred_multidomain_dirs.R](create_plots_from_tarred_multidomain_dirs.R) can make basic png images of various flow maxima (stage, speed, flux) from the tarred multidomain directories.
* [make_rasters.R](make_rasters.R) is useful for making rasters from a single model run.
* [load_balance_script.R](load_balance_script.R) can be run from inside a multidomain directory, and will produce a file `load_balance_file.txt` which tries to distribute domains among MPI images to equidistribute the work from that run. 
    * The files in [./load_balance_files](./load_balance_files) were produced this way (but apply to different model setups and core-counts). 
* [make_domains_shapefile.R](make_domains_shapefile.R) can make a shapefile depicting the multidomain layout. 
* [make_initial_conditions_complex_historical_events.R](make_initial_conditions_complex_historical_events.R) makes an input file for the Sumatra 2004 source-inversion of [Fujii and Satake, 2007](https://doi.org/10.1785/0120050613), which has a different time-history of rupture for each unit-source.
