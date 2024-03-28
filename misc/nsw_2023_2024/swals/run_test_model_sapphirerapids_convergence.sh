#!/bin/bash
#PBS -P w85
#PBS -q normalsr
#PBS -l walltime=00:30:00
#PBS -lmem=6000GB
#PBS -lncpus=1248
#PBS -l wd
#PBS -l storage=scratch/w85+gdata/w85

#source SWALS_ifort_modules_2023_llvm_B.sh
source SWALS_ifort_modules_2023_B.sh

#OMP_NUM_THREADS=13 mpiexec -np 96 --map-by ppr:4:socket:PE=13 ./model_sapphirerapids \
mpiexec -np 96 --map-by numa:SPAN --bind-to numa -x OMP_NUM_THREADS=13 -x OMP_PROC_BIND=true ./model_sapphirerapids \
    '' multidomain_design_control_NNL4_1arcminoffshore_defaultloadbalance_CONVERGENCE.nml \
    stationary_test_with_LHI_NNL4_CONVERGENCE test_load_balance 1.0 > outfile.log
