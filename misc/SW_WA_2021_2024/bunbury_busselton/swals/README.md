Tsunami modelling codes for the Bunbury/Busselton region
--------------------------------------------------------

This folder contains code to run the tsunami models, including:
* Scenarios like historical events [(see here)](../sources/like_historic/)
* Random PTHA18 scenarios [(see here)](../sources/hazard/)

The key codes are described below; consult comments in the code itself for further documentation.

## Three versions of this model

1. The tsunami model was initially setup to assume the Bunbury floodgate was open. Notice the Bunbury floodgate is created in the folder [../breakwalls](../breakwalls), but it isn't added to the default list of breakwalls. Instead there is a parameter `logical, parameter :: close_bunbury_floodgate = .FALSE.` in [model_local_routines.f90](model_local_routines.f90) that determines if it is used.

2. Later (mid 2023) the model was re-run with the Bunbury floodgate closed, by setting the parameter `close_bunbury_floodgate=.true.`.

3. Later (April 2024) the model was re-run with improved data near the Vasse diversion drain (breakwalls and elevation patches). A new parameter `logical, parameter :: include_elevation_updates_2024 = .TRUE.` was added to the code so this could be done in a backward compatible manner. For this model we had the Bunbury floodgate open.

## Hydrodynamic model setup and compilation

* [model.f90](model.f90) and [model_local_routines.f90](model_local_routines.f90) are used to setup the hydrodynamic model.
    * [make_model_ifort](make_model_ifort) was used to compile the model on NCI, after loading the modules in [SWALS_ifort_modules_2022.sh](SWALS_ifort_modules_2022.sh) 
        * `source SWALS_ifort_modules_2022.sh; make -B -f make_model_ifort`
    * The model requires information on the multidomain layout ([created here](../multidomain_design/)) and the initial conditions. The latter include [the initial water-surface perturbation](../sources/) and [the elevation grids](../elevation/) and [the breakwall geometry](../breakwalls).

## Creation of qsub scripts to run the hydrodynamic model for random PTHA scenarios

* [create_random_ptha_qsub_scripts_sealevel60cm.R](create_random_ptha_qsub_scripts_sealevel60cm.R) is used to make qsub scripts which run the hydrodynamic model for all the random scenarios (for the full-resolution runs), and also `tar` the resulting output folders (to prevent creation of too many files on NCI). 
    * The script can be run with 
        * `Rscript create_random_ptha_qsub_scripts_sealevel60cm.R`.
    * It produces approximately 50 qsub scripts which can be separately submitted. One can also submit all the scripts at once using shell commands such as: 
        * `for i in run_BunburyBusselton_sealevel60cm_*.sh; do qsub $i; mv $i submitted_qsub_jobs; done`. 

* [DEBUG_run_with_old_nesting.sh](DEBUG_run_with_old_nesting.sh) was used to run a scenario that went unstable after I updated to the new Busselton data in March 2023. Of the 369 scenarios, only one was unstable. As a workaround I compiled SWALS with the old nesting algorithm using [make_model_ifort_OLD_NESTING](make_model_ifort_OLD_NESTING), and reran the scenario. It did not go unstable with the old nesting approach (in general, that approach isn't better -- but nesting instabilities are rare and fragile, so often changing the details will fix issues). I manually removed the old run, copied the `multidomain_*0001.log` file to the right place, and re-ran the script to create rasters. At this point the model appears just like the other scenarios; for our purposes, no more workarounds are required.

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
* There are some scripts that are designed to copy [gauges](make_folders_and_copy_gauges.R) or [rasters](make_folders_and_copy_rasters.R) from NCI to a local machine. They need modification to point to a file containing your NCI username and password.
* [set_max_stage_to_NA_where_max_flux_is_zero.R](set_max_stage_to_NA_where_max_flux_is_zero.R) modifies max-stage rasters to set the stage to zero when the max-flux is zero. This is a cosmetic workaround for sites in the zone of initial earthquake subsidence (i.e. far from the site where we model inundation). This was done for a couple of graphics because the SWALS max-stage calculation considers the max-stage over time. If a time-varying source is used, then in areas with subsidence, it is possible for the max-stage to simply reflect the pre-subsidence topography even if the site is dry. The SWALS elevation output gives elevation at the last timestep, and so our native processing could spuriously make some "dry" subsiding sites appear wet. A simple workaround is to zero the max-stage in sites where the max-flux is zero. A better workaround (future) might be to adjust the SWALS max-stage calculation to ignore dry cells, if that doesn't have other negative consequences.
