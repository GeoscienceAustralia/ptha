#!/bin/bash
#PBS -P w85
#PBS -q normal
#PBS -l walltime=15:00:00
#PBS -lmem=2280GB
#PBS -lncpus=576
#PBS -l wd
#PBS -l storage=scratch/w85+gdata/w85

source SWALS_ifort_modules_2021.sh

# Regular case
OMP_NUM_THREADS=12 mpiexec -np 48 --map-by ppr:2:socket:PE=12 ./model_meshrefine8 '../sources/like_historic/FujiSatake2007/Fuji_andaman2004_unit_sources_SUM_KAJIURA_SMOOTHED.tif' Fuji_andaman2004_24hrs_Mandurah2Geraldton_meshrefine8 full '../multidomain_design/domains_0.5_0.125/first_level_nesting.csv' '../multidomain_design/domains_0.5_0.125/second_level_nesting.csv' 'load_balance_files/load_balance_011121_9x_0.5_0.125_48MPI_meshrefine8.txt' 0.0 > outfile.log


# Test load balance
#OMP_NUM_THREADS=12 mpiexec -np 48 --map-by ppr:2:socket:PE=12 ./model_meshrefine8 '../sources/like_historic/FujiSatake2007/Fuji_andaman2004_unit_sources_SUM_KAJIURA_SMOOTHED.tif' Fuji_andaman2004_24hrs_Mandurah2Geraldton_meshrefine8 test_load_balance '../multidomain_design/domains_0.5_0.125/first_level_nesting.csv' '../multidomain_design/domains_0.5_0.125/second_level_nesting.csv' '../multidomain_design/domains_0.5_0.125/load_balance_default.txt' 0.0 > outfile.log
