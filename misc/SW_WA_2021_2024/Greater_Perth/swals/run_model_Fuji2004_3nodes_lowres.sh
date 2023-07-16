#!/bin/bash
#PBS -P w85
#PBS -q normal
#PBS -l walltime=02:00:00
#PBS -lmem=570GB
#PBS -lncpus=144
#PBS -l wd
#PBS -l storage=scratch/w85+gdata/w85

source SWALS_ifort_modules_2022.sh

#
# Case using a third level of nesting
#

# New model, low-res, testing original load balancing
#OMP_NUM_THREADS=6 mpiexec -np 24 --map-by ppr:4:socket:PE=6 ./model '../sources/like_historic/FujiSatake2007/Fuji_andaman2004_unit_sources_SUM_KAJIURA_SMOOTHED.tif' Fuji_andaman2004_24hrs_domain010322_lowres test_load_balance '../multidomain_design/domains_010322_0.5_0.166666666666667_0.0333333333333333/first_level_nesting_edited.csv' '../multidomain_design/domains_010322_0.5_0.166666666666667_0.0333333333333333/second_level_nesting_edited.csv' '../multidomain_design/domains_010322_0.5_0.166666666666667_0.0333333333333333/third_level_nesting_edited.csv' 'load_balance_files/load_balance_020322_0.5_0.166666666666667_0.0333333333333333_24MPI.txt' 0.0 > outfile.log

## New model, testing custom balancing
# OMP_NUM_THREADS=6 mpiexec -np 24 --map-by ppr:4:socket:PE=6 ./model '../sources/like_historic/FujiSatake2007/Fuji_andaman2004_unit_sources_SUM_KAJIURA_SMOOTHED.tif' Fuji_andaman2004_24hrs_domain010322_lowres test_load_balance '../multidomain_design/domains_010322_0.5_0.166666666666667_0.0333333333333333/first_level_nesting_edited.csv' '../multidomain_design/domains_010322_0.5_0.166666666666667_0.0333333333333333/second_level_nesting_edited.csv' '../multidomain_design/domains_010322_0.5_0.166666666666667_0.0333333333333333/third_level_nesting_edited.csv' 'load_balance_files/load_balance_120322_0.166666666666667_0.0333333333333333_24MPI_lowres.txt' 0.0 > outfile.log

## New model, serious
#OMP_NUM_THREADS=6 mpiexec -np 24 --map-by ppr:4:socket:PE=6 ./model '../sources/like_historic/FujiSatake2007/Fuji_andaman2004_unit_sources_SUM_KAJIURA_SMOOTHED.tif' Fuji_andaman2004_24hrs_domain010322_lowres full '../multidomain_design/domains_010322_0.5_0.166666666666667_0.0333333333333333/first_level_nesting_edited.csv' '../multidomain_design/domains_010322_0.5_0.166666666666667_0.0333333333333333/second_level_nesting_edited.csv' '../multidomain_design/domains_010322_0.5_0.166666666666667_0.0333333333333333/third_level_nesting_edited.csv' 'load_balance_files/load_balance_120322_0.166666666666667_0.0333333333333333_24MPI_lowres.txt' 0.0 > outfile.log

## New model, time-varying forcing
OMP_NUM_THREADS=6 mpiexec -np 24 --map-by ppr:4:socket:PE=6 ./model 'FujiSatake2007_time_varying_forcing_realistic.csv' Fuji_andaman2004_24hrs_domain010322_lowres_timevaryingRealistic full '../multidomain_design/domains_010322_0.5_0.166666666666667_0.0333333333333333/first_level_nesting_edited.csv' '../multidomain_design/domains_010322_0.5_0.166666666666667_0.0333333333333333/second_level_nesting_edited.csv' '../multidomain_design/domains_010322_0.5_0.166666666666667_0.0333333333333333/third_level_nesting_edited.csv' 'load_balance_files/load_balance_120322_0.166666666666667_0.0333333333333333_24MPI_lowres.txt' 0.0 > outfile.log

