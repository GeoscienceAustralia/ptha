#!/bin/bash
#PBS -P w85
#PBS -q normalsr
#PBS -l walltime=00:20:00
#PBS -lmem=1000GB
#PBS -lncpus=208
#PBS -l wd
#PBS -l storage=scratch/w85+gdata/w85

#source SWALS_ifort_modules_2023_llvm_B.sh
source SWALS_ifort_modules_2023_B.sh

mpiexec -np 16 --map-by numa:SPAN --bind-to numa -x OMP_NUM_THREADS=13 -x OMP_PROC_BIND=true ./model_sapphirerapids \
    '' \
    multidomain_design_control_NNL4_1arcminoffshore.nml \
    stationary_test_with_LHI_NNL4 test_load_balance 1.0 > outfile.log
