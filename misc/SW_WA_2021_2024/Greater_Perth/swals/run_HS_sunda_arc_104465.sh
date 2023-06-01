#!/bin/bash
#PBS -P w85
#PBS -q normal
#PBS -l walltime=5:00:00
#PBS -lmem=1140GB
#PBS -lncpus=288
#PBS -l wd
#PBS -l storage=scratch/w85+gdata/w85

source SWALS_ifort_modules_2021.sh

OMP_NUM_THREADS=12 mpiexec -np 24 --map-by ppr:2:socket:PE=12 ./model '../../../australia_wide/clean_version/ptha18_scenarios_random/set_range_of_mw_and_centroid/sumatra2004/sumatra2004_heterogeneous_slip_104465_count_1.tif' HS_sunda_arc_104465_Mandurah2Geraldton full '../multidomain_design/domains_0.5_0.125/first_level_nesting.csv' '../multidomain_design/domains_0.5_0.125/second_level_nesting.csv' 'load_balance_files/load_balance_120821_9x_0.5_0.125_24MPI.txt' 0.0 > outfile.log


