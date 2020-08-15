# This folder contains code to run various SWALS tsunami models for the inversions. 
---------------------------------------------------------------------------------

## The key files are:

* `model.f90` and `model_local_routines.f90`  -- Tsunami model setup codes. Designed for this study only. See comments therein for documentation.

* `point_gauges_combined.csv` -- Coordinates of all sites at which the tsunami model stores gauge time-series. A small subset of these have observational data and are used in further processing; many others are points from the 2018 PTHA.

* `make_model_ifort`, `SWALS_ifort_modules.sh` -- These are used to load software and compile `model.f90` on NCI's Gadi machine. They are run with `source SWALS_ifort_modules.sh; make -B -f make_model_ifort`.

* Files like `load_balance_files/load_balance_XXXXX.txt` -- These files tell SWALS how to split the domains for distributed memory parallel runs, in a way that is efficient with our machine, software and compiler options. Further information on the load balance files is presented below.

* Files like `run_*.sh` -- Submit jobs to the PBS queuing system. 

* Code in `./plots` -- Compare observations and models at gauges. Importantly these scripts also save the model and data to an RDS file, which is useful for subsequent analysis. We run `./plots/plot_all.R` after all the SWALS jobs have completed to make the plots for all jobs..

* `copy_gauges_to_analysis_directory.R` -- After all the SWALS models have been run, AND the codes in `./plots` have been run, then run this script to copy the gauge and relevant log data to a folder in `../analysis` for subsequent processing.

* Other misc. R scripts are present but not essential.

## Procedure to run the tsunami models

Assuming all the dependent files are set-up correctly (i.e. elevation data is under the `../elevation` folder, breakwalls are in the `../breakwalls` folder, sources are in the `../sources` folder, etc), then on NCI's gadi machine the models can be run like this:

* Load the required compilers / libraries etc:

  `source SWALS_ifort_modules.sh`

* Compile the code

    `make -B -f make_model_ifort`

* At this point the code is ready to run. Below are the scripts we used for this study, which use the PBS queueing system on the Gadi supercomputer. These would need to be edited to run on a different machine.

```

    # To check the effect of the rise-time (offshore only models, one inversion only)
    # Mathematically for linear equations, a rise-time is equivalent to smoothing in time 
    # (and this is tested in the SWALS validation suite).
    qsub run_inversions_tohoku2011_2nodes_offshoreOnly_rise_time_sensitivity.sh

    # To support a convergence test -- the job is like the above rise_time=0 case, with 2x higher resolution.
    qsub run_inversions_tohoku2011_4nodes_offshoreOnly_HIGHER_RES.sh

    # Run "purely offshore' models for a range of scenarios. I used this for preliminary study of 
    # energy decay rates (noting the global & DART behaviour is very similar to runs that have nested 
    #regions around Australia).
    qsub run_inversions_with_DART_obs_linear_only_4nodes.sh
    qsub run_inversions_with_DART_obs_linear_with_manning_4nodes.sh

    # This runs 'test-jobs' that can be used to make load-balance files (the ones used are in ./load_balance_files/)
    qsub run_inversions_tohoku2011_8nodes_TEST_RUNS_FOR_LOAD_BALANCING.sh
    #
    # When the above 4 jobs have completed, one can create a load-balance file for each combination of [NSW OR australia] and
    # [linear_with_manning OR linear_with_linear_friction]. Note load-balance files for linear_with_linear_friction can also be used
    # for linear_with_no_friction. The procedure to make the files is
    #    
    #     1) On NCI, get the R modules
    #        source R_361_NCI_modules.sh 
    #     2) For each model change directory into its multidomain directory. For instance, the 'linear_with_manning + NSW' case is like:
    #        cd OUTPUTS/Tohoku2011_Romano2015-risetime_0.0-test_load_balance-linear_with_manning-0.035-highres_NSW/RUN_20200604_180848234/
    #     3) Start R
    #        # Source the SWALS/plot.R scripts -- for instance on NCI I do this from inside R
    #        source('/g/data/w85/tsunami/CODE/gadi/ptha/propagation/SWALS/plot.R')
    #        # Create the load balance file. This will print information on the benefits expected from load balancing
    #        make_load_balance_partition('.')
    #        # If you are satisfied with the expected benefit, then copy the resulting load_balance_partiton.txt file to
    #        # the directory load_balance_files (the latter is inside the directory holding this README.md). Give it a name that
    #        # matches the job-type, e.g. for the case with 'linear_with_manning + NSW', do
    #        file.copy('load_balance_partition.txt', '../../../load_balance_files/load_balance_manning_offshore_NSW_8nodes_32ranks.txt')
    # Repeat the above for the other cases
    # Then we can run the main jobs
   
    # Run all cases that have no-friction in the offshore solver
    qsub run_inversions_full_linear_with_no_friction_8nodes.sh

    # Run all cases that have manning-friction in the offshore solver
    qsub run_inversions_full_linear_with_manning_8nodes.sh

    # Run all cases that have a constant linear friction (1e-05) in the offshore solver
    qsub run_inversions_full_linear_with_linear_friction_8nodes.sh
    
    # Run all cases that have a delayed constant linear friction (1e-05) in the offshore solver
    qsub run_inversions_full_linear_with_delayed_linear_friction_8nodes.sh

    # Run all cases that have a reduced constant linear friction (1/(36*3600)) in the offshore solver
    qsub run_inversions_full_linear_with_reduced_linear_friction_8nodes.sh

    # Run one of the models with a refined mesh, to check convergence
    qsub run_yamakazi_HIGHRES_30nodes.sh

    # Run another refined mesh code to check convergence
    run_inversions_Fujii07_HIGHRES_30nodes.sh

    # Run a model using leapfrog-nonlinear offshore, for comparison with "linear+manning" offshore
    run_yamakazi_leapfrognonlinearoffshore_8nodes.sh
```

* At this point run the code in `./plots` by changing into that directory and running `Rscript plot_all.R`. This will create a lot of plots and also RDS files with observations and models at gauges.

* Next make png files with the max stage and elevation for each model with `Rscript plot_max_stage_and_elevation.R` (note the hard-coded assumption of 24 cores).

* To copy the gauge results and plots to a folder in ../analysis for further processing, do `Rscript copy_gauges_to_analysis_directory.R`

## Further information on the load balance file format

The `model.f90` code creates a multidomain object `md` that contains an array of domains `md%domains(:)`. Each domain corresponds to a structured grid where we solve the shallow water equations. To run in distributed-memory parallel, the code splits each domain up into a certain number of pieces and distributes them among MPI ranks. Ideally we will partition the domains in a way that minimises the total runtime (so roughly equal work is given to each MPI process). While SWALS can create a default partition, this may be far from optimal. Herein we specify the partition to improve the code efficiency, using files in `./load_balance_files/load_balance_XXX.txt` files. 

Here we explain the format of the load balance files ([see this example](./load_balance_files/load_balance_linear_offshore_NSW_8nodes_32ranks.txt)). Each line corresponds to a domain, and contains several integers. For example if there are 16 domains in the multidomain (`size(md%domains) == 16`) then the file must have 16 lines. On each line, every integer corresponds to a coarray-image (equal to the "mpi-rank + 1"). SWALS will split that domain among the specified coarray-images in roughly equal pieces. For example; if the 8th row contains 4 integers then md%domains(8) will be split into 4 roughly equal pieces among the specified coarray images. Integers may be repeated, to assign more than one piece of the domain to an image.

How can we create a good load-balance file? In general the best partition will depend on your machine, compiler, compiler options, and perhaps details of the model input data. Let's assume all these are fixed. The first step is to manually create a default load-balance file by guessing how many pieces each domain should be split into. Assign nominal integer values to the images, ranging from 1 up to the number of images you will use to run the model ([here is an example](./load_balance_files/load_balance_default_NSW_8nodes_32mpi.txt)). The default load-balance file is used to run a short simulation (passing the `test_load_balance` option to `./model`). During this run, the time required for each domain is recorded in the model log files. The log files can be used by the SWALS function `make_load_balance_partition` (this is in the ptha repository at [ptha/propagation/SWALS/plot.R](https://github.com/GeoscienceAustralia/ptha/blob/master/propagation/SWALS/plot.R)) to a make an improved load-balance file, which attempts to minimise the runtime using a Greedy approach to load balancing. 

`make_load_balance_partition` will print an estimate of how much time will be saved and how much load-imbalance will remain. If the balance is good, we can use the resulting `load_balance_partition.txt` file for real applications (that's how [this_example](./load_balance_files/load_balance_linear_offshore_NSW_8nodes_32ranks.txt) was made). Depending on our original guess as to how many pieces each domain should be partitioned into, it is possible we cannot make a good load-balance file. For instance, this might happen if one or a few domain pieces dominate the computational effort. In those cases we can revise our guess in the default load balance file and try again.
