#!/bin/bash
#PBS -P w85
#PBS -q normal
#PBS -l walltime=10:00:00
#PBS -lmem=384GB
#PBS -lncpus=96
#PBS -l wd
#PBS -l storage=scratch/n74

source SWALS_ifort_modules.sh
ulimit -s unlimited



# Tohoku, PTHA18 HS scenario 46994 
OMP_NUM_THREADS=24 OMP_PROC_BIND=true mpiexec -np 4 --map-by ppr:1:socket:PE=24 ./model ../sources/like_historic/kurilsjapan_46994_stochastic_gauge_summary_stats_session_kurilsjapan_tohoku_2011_03_11_Mw9.1.Rdata.tif 0 0.0 Tohoku2011_validation_PTHA18_HS_46994 'full' '' linear_with_manning 0.035 tonga pacific few_grids > outfile.log
