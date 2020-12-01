#!/bin/bash
#PBS -P w85
#PBS -q normal
#PBS -l walltime=04:00:00
#PBS -lmem=384GB
#PBS -lncpus=96
#PBS -l wd
#PBS -l storage=scratch/n74

source SWALS_ifort_modules.sh
ulimit -s unlimited


# Tonga 2006 scenario
OMP_NUM_THREADS=24 OMP_PROC_BIND=true mpiexec -np 4 --map-by ppr:1:socket:PE=24 ./model ../sources/like_historic/kermadectonga2_26849_variable_uniform_gauge_summary_stats_session_kermadectonga2_tonga_2006_05_03_Mw8.0.Rdata.tif 0 0.0 Tonga2006_validation_PTHA18_VAUS_26849_load_balance 'full' 'load_balance_partition.txt' linear_with_manning 0.035 tonga regional few_grids > outfile.log

