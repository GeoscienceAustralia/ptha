#!/bin/bash
#PBS -P w85
#PBS -q normalsr
#PBS -l walltime=08:00:00
#PBS -lmem=1000GB
#PBS -lncpus=208
#PBS -l wd
#PBS -l storage=scratch/w85+gdata/w85

#source SWALS_ifort_modules_2023_llvm_B.sh
source SWALS_ifort_modules_2023_B.sh

#OMP_NUM_THREADS=13 mpiexec -np 16 --map-by ppr:4:socket:PE=13 ./model_sapphirerapids \
mpiexec -np 16 --map-by numa:SPAN --bind-to numa -x OMP_NUM_THREADS=13 -x OMP_PROC_BIND=true ./model_sapphirerapids \
    '../sources/like_historic/Puysegur2009/PTHA18_HS_1567/puysegur2_1567_stochastic_gauge_summary_stats_session_puysegur2_puysegur_2009_07_15_Mw7.8.Rdata.tif' \
    'multidomain_design_control_NNL4_1arcminoffshore.nml' \
    run_Puysegur2009_PTHA18HS1567_1arcminoffshore full 0.0 > outfile.log
