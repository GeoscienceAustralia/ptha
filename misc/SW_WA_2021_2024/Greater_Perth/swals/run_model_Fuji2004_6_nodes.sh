#!/bin/bash
#PBS -P w85
#PBS -q normal
#PBS -l walltime=04:00:00
#PBS -lmem=1140GB
#PBS -lncpus=288
#PBS -l wd
#PBS -l storage=scratch/w85+gdata/w85

source SWALS_ifort_modules_2022.sh

#
# Case using a third level of nesting
#

## New model, testing default balancing
#OMP_NUM_THREADS=12 mpiexec -np 24 --map-by ppr:2:socket:PE=12 ./model '../sources/like_historic/FujiSatake2007/Fuji_andaman2004_unit_sources_SUM_KAJIURA_SMOOTHED.tif' Fuji_andaman2004_24hrs_domain181221 test_load_balance '../multidomain_design/domains_181221_0.5_0.166666666666667_0.0333333333333333/first_level_nesting_edited.csv' '../multidomain_design/domains_181221_0.5_0.166666666666667_0.0333333333333333/second_level_nesting_edited.csv' '../multidomain_design/domains_181221_0.5_0.166666666666667_0.0333333333333333/third_level_nesting_edited.csv' '../multidomain_design/domains_181221_0.5_0.166666666666667_0.0333333333333333/load_balance_default.txt' 0.0 > outfile.log

## New model, testing custom balancing
#OMP_NUM_THREADS=12 mpiexec -np 24 --map-by ppr:2:socket:PE=12 ./model '../sources/like_historic/FujiSatake2007/Fuji_andaman2004_unit_sources_SUM_KAJIURA_SMOOTHED.tif' Fuji_andaman2004_24hrs_domain181221 test_load_balance '../multidomain_design/domains_181221_0.5_0.166666666666667_0.0333333333333333/first_level_nesting_edited.csv' '../multidomain_design/domains_181221_0.5_0.166666666666667_0.0333333333333333/second_level_nesting_edited.csv' '../multidomain_design/domains_181221_0.5_0.166666666666667_0.0333333333333333/third_level_nesting_edited.csv' 'load_balance_files/load_balance_181221_0.5_0.166666666666667_0.0333333333333333_24MPI.txt' 0.0 > outfile.log

## New model, serious
#OMP_NUM_THREADS=12 mpiexec -np 24 --map-by ppr:2:socket:PE=12 ./model '../sources/like_historic/FujiSatake2007/Fuji_andaman2004_unit_sources_SUM_KAJIURA_SMOOTHED.tif' Fuji_andaman2004_24hrs_domain181221_newProcessDataToSend full '../multidomain_design/domains_181221_0.5_0.166666666666667_0.0333333333333333/first_level_nesting_edited.csv' '../multidomain_design/domains_181221_0.5_0.166666666666667_0.0333333333333333/second_level_nesting_edited.csv' '../multidomain_design/domains_181221_0.5_0.166666666666667_0.0333333333333333/third_level_nesting_edited.csv' 'load_balance_files/load_balance_181221_0.5_0.166666666666667_0.0333333333333333_24MPI.txt' 0.0 > outfile.log

## New model, time-varying forcing
OMP_NUM_THREADS=12 mpiexec -np 24 --map-by ppr:2:socket:PE=12 ./model 'FujiSatake2007_time_varying_forcing_realistic.csv' Fuji_andaman2004_24hrs_domain010322_timevaryingRealistic full '../multidomain_design/domains_010322_0.5_0.166666666666667_0.0333333333333333/first_level_nesting_edited.csv' '../multidomain_design/domains_010322_0.5_0.166666666666667_0.0333333333333333/second_level_nesting_edited.csv' '../multidomain_design/domains_010322_0.5_0.166666666666667_0.0333333333333333/third_level_nesting_edited.csv' 'load_balance_files/load_balance_020322_0.5_0.166666666666667_0.0333333333333333_24MPI.txt' 0.0 > outfile.log

## New model, instantaneous version of time-varying forcing for testing
#OMP_NUM_THREADS=12 mpiexec -np 24 --map-by ppr:2:socket:PE=12 ./model 'FujiSatake2007_time_varying_forcing_instantaneous.csv' Fuji_andaman2004_24hrs_domain181221_timevaryingInstantaneous full '../multidomain_design/domains_181221_0.5_0.166666666666667_0.0333333333333333/first_level_nesting_edited.csv' '../multidomain_design/domains_181221_0.5_0.166666666666667_0.0333333333333333/second_level_nesting_edited.csv' '../multidomain_design/domains_181221_0.5_0.166666666666667_0.0333333333333333/third_level_nesting_edited.csv' 'load_balance_files/load_balance_181221_0.5_0.166666666666667_0.0333333333333333_24MPI.txt' 0.0 > outfile.log
