#!/bin/bash
#PBS -P w85
#PBS -q normalsr
#PBS -l walltime=15:00:00
#PBS -lmem=6000GB
#PBS -lncpus=1248
#PBS -l wd
#PBS -l storage=scratch/w85+gdata/w85

#source SWALS_ifort_modules_2023_llvm_B.sh
source SWALS_ifort_modules_2023_B.sh

#OMP_NUM_THREADS=13 mpiexec -np 96 --map-by ppr:4:socket:PE=13 ./model_sapphirerapids \
mpiexec -np 96 --map-by numa:SPAN --bind-to numa -x OMP_NUM_THREADS=13 -x OMP_PROC_BIND=true ./model_sapphirerapids \
    '../../../../east_australian_coast_2021_02/sources/kt_multi_site_selection/kermadectonga2_stochastic_43731_9.4_10000.tif' \
    'multidomain_design_control_NNL4_1arcminoffshore_12hrs_final_CONVERGENCE.nml' \
    run_kt43731_12h_final_NNL4_CONVERGENCE full 1.1 > outfile.log
