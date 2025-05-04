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

mpiexec -np 16 --map-by numa:SPAN --bind-to numa -x OMP_NUM_THREADS=13 -x OMP_PROC_BIND=true ./model_sapphirerapids \
    '../sources/like_historic/Chile2014/southamerica_80427_stochastic_varyMu_gauge_summary_stats_session_southamerica_southamerica_2014_04_01_Mw8.2_varyMu.Rdata.tif' \
    'multidomain_design_control_NNL4_1arcminoffshore.nml' \
    run_Chile2014_PTHA18_HS80427_1arcminoffshore full 0.0 > outfile.log
