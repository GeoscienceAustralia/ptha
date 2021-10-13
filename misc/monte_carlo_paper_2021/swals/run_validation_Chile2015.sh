#!/bin/bash
#PBS -P w85
#PBS -q normal
#PBS -l walltime=10:00:00
#PBS -lmem=384GB
#PBS -lncpus=96
#PBS -l wd
#PBS -l storage=gdata/n74+gdata/w85

source SWALS_ifort_modules_2021.sh
ulimit -s unlimited

# Chile 2015, PTHA18 scenario 101114
OMP_NUM_THREADS=24 OMP_PROC_BIND=true mpiexec -np 4 --map-by ppr:1:socket:PE=24 ./model '../sources/like_historic/southamerica_101114_stochastic_gauge_summary_stats_session_southamerica_southamerica_2015_09_16_Mw8.3.Rdata.tif' 0 0.0 Chile2015_validation_PTHA18_HS_101114_source 'full' '' linear_with_manning 0.035 tonga pacific few_grids > outfile.log


