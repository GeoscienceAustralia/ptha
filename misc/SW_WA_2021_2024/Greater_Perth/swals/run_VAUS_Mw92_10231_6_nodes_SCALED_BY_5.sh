#!/bin/bash
#PBS -P w85
#PBS -q normal
#PBS -l walltime=04:00:00
#PBS -lmem=1140GB
#PBS -lncpus=288
#PBS -l wd
#PBS -l storage=scratch/w85+gdata/w85

source SWALS_ifort_modules_2021.sh

OMP_NUM_THREADS=12 mpiexec -np 24 --map-by ppr:2:socket:PE=12 ./model '../sources/testing/VAUS_Mw92_sunda_arc/Mw92_scenarios_similar_size/Mw92_initial_condition_10231.tif' VAUS_10231_SCALEDBY5_domains010322_newProcessDataToSend full '../multidomain_design/domains_010322_0.5_0.166666666666667_0.0333333333333333/first_level_nesting_edited.csv' '../multidomain_design/domains_010322_0.5_0.166666666666667_0.0333333333333333/second_level_nesting_edited.csv' '../multidomain_design/domains_010322_0.5_0.166666666666667_0.0333333333333333/third_level_nesting_edited.csv' 'load_balance_files/load_balance_020322_0.5_0.166666666666667_0.0333333333333333_24MPI.txt' 0.6 > outfile.log

