#!/bin/bash
#PBS -P w85
#PBS -q normalsr
#PBS -l walltime=10:00:00
#PBS -lmem=1000GB
#PBS -lncpus=208
#PBS -l wd
#PBS -l storage=scratch/w85+gdata/w85

#source SWALS_ifort_modules_2023_llvm_B.sh
source SWALS_ifort_modules_2023_B.sh

#OMP_NUM_THREADS=13 mpiexec -np 16 --map-by ppr:4:socket:PE=13 ./model_sapphirerapids \
mpiexec -np 16 --map-by numa:SPAN --bind-to numa -x OMP_NUM_THREADS=13 -x OMP_PROC_BIND=true ./model_sapphirerapids \
    '../sources/test/small_source_solomon2007_divided_by_1000/small_scenario_solomon2007_divided_by_1000.tif' \
    'multidomain_design_control_NNL4_1arcminoffshore_timegrids.nml' \
    run_small_source_1arcminoffshore full 0.0 > outfile.log