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



# Chile 2015, inversion based on Williamson et al 
OMP_NUM_THREADS=24 OMP_PROC_BIND=true mpiexec -np 4 --map-by ppr:1:socket:PE=24 ./model ../sources/like_historic/Williamson_chile_sources_SUM_KAJIURA_SMOOTHED.tif 0 0.0 Chile2015_validation_Williamson_source_inversion 'full' '' linear_with_manning 0.035 tonga pacific few_grids > outfile.log
