#!/bin/bash
#PBS -P w85
#PBS -q normal
#PBS -l walltime=04:00:00
#PBS -lmem=384GB
#PBS -lncpus=96
#PBS -l wd
#PBS -l storage=gdata/n74+gdata/w85

source SWALS_ifort_modules_2021.sh
ulimit -s unlimited

# Tonga 2009 scenario
#OMP_NUM_THREADS=24 OMP_PROC_BIND=true mpiexec -np 4 --map-by ppr:1:socket:PE=24 ./model ../sources/like_historic/kermadectonga2_24563_stochastic_gauge_summary_stats_session_kermadectonga2_samoa_2009_09_29_Mw8.1.Rdata.tif 0 0.0 Tonga2009_validation_PTHA18_HS_24563_load_balance 'full' 'load_balance_partition.txt' linear_with_manning 0.035 tonga regional few_grids > outfile.log

OMP_NUM_THREADS=24 OMP_PROC_BIND=true mpiexec -np 4 --map-by ppr:1:socket:PE=24 ./model ../sources/like_historic/kermadectonga2_28606_stochastic_gauge_summary_stats_session_kermadectonga2_samoa_2009_09_29_Mw8.1.Rdata.tif 0 0.0 Tonga2009_validation_PTHA18_HS_28606_load_balance 'full' 'load_balance_partition.txt' linear_with_manning 0.035 tonga regional few_grids > outfile.log

#OMP_NUM_THREADS=24 OMP_PROC_BIND=true mpiexec -np 4 --map-by ppr:1:socket:PE=24 ./model ../sources/like_historic/kermadectonga2_28678_stochastic_gauge_summary_stats_session_kermadectonga2_samoa_2009_09_29_Mw8.1.Rdata.tif 0 0.0 Tonga2009_validation_PTHA18_HS_28678_load_balance 'full' 'load_balance_partition.txt' linear_with_manning 0.035 tonga regional few_grids > outfile.log
