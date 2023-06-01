#!/bin/bash
#PBS -P w85
#PBS -q normal
#PBS -l walltime=06:00:00
#PBS -lmem=1140GB
#PBS -lncpus=288
#PBS -l wd
#PBS -l storage=scratch/w85+gdata/w85

source SWALS_ifort_modules_2021.sh

#
# Case using a third level of nesting
#

OMP_NUM_THREADS=12 mpiexec -np 24 --map-by ppr:2:socket:PE=12 ./model '../sources/like_historic/Fujii2021/Fuji21_andaman2004_unit_sources_SUM_KAJIURA_SMOOTHED.tif' Fujii2021_andaman2004_24hrs_threeLevelNesting_ThreeMinRiseTime full '../multidomain_design/domains_0.5_0.166666666666667_0.0333333333333333/first_level_nesting.csv' '../multidomain_design/domains_0.5_0.166666666666667_0.0333333333333333/second_level_nesting.csv' '../multidomain_design/domains_0.5_0.166666666666667_0.0333333333333333/third_level_nesting.csv' 'load_balance_files/load_balance_181121_0.5_0.166666666666667_0.0333333333333333_24MPI.txt' 0.0 > outfile.log
