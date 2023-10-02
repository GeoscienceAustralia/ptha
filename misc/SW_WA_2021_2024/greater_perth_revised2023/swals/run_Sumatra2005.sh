#!/bin/bash
#PBS -P w85
#PBS -q normalsr
#PBS -l walltime=04:00:00
#PBS -lmem=1000GB
#PBS -lncpus=208
#PBS -l wd
#PBS -l storage=scratch/w85+gdata/w85

source SWALS_ifort_modules_2023_llvm.sh

# Run with sea level = 0.0 m
OMP_NUM_THREADS=13 mpiexec -np 16 --map-by ppr:4:socket:PE=13 ./model_sapphirerapids \
    ../sources/like_historic/sumatra2005/Fuji_nias2005_unit_sources_SUM_KAJIURA_SMOOTHED.tif \
    multidomain_design_control_NNL4_defaultres.nml \
    Sumatra2005_FujiiEtAl2020 \
    full \
    0.0 > outfile.log
