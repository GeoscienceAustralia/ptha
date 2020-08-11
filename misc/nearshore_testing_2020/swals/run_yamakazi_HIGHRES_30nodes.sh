#!/bin/bash
#PBS -P w85
#PBS -q normal
#PBS -l walltime=05:00:00
#PBS -lmem=5700GB
#PBS -lncpus=1440
#PBS -l wd
#PBS -l storage=scratch/w85+gdata/w85

source SWALS_ifort_modules.sh
ulimit -s unlimited

# This is for convergence testing -- compiled with "mesh_refine=8_ip"

## For load balance file creation
#OMP_NUM_THREADS=12 OMP_PROC_BIND=true mpiexec -np 120 --map-by ppr:2:socket:PE=12 ./model ../sources/Tohoku2011/YamakaziEtAl2018/yamakazi18_Tohoku11_sources_SUM_KAJIURA_SMOOTHED.tif 0 Tohoku2011_YamakaziEtAl2018_HIGHRES test_load_balance load_balance_files/load_balance_default_NSW_30nodes_120mpi.txt linear_with_manning 0.035 NSW > outfile.log

# Full
OMP_NUM_THREADS=12 OMP_PROC_BIND=true mpiexec -np 120 --map-by ppr:2:socket:PE=12 ./model ../sources/Tohoku2011/YamakaziEtAl2018/yamakazi18_Tohoku11_sources_SUM_KAJIURA_SMOOTHED.tif 0 Tohoku2011_YamakaziEtAl2018_HIGHRES full load_balance_files/load_balance_manning_offshore_NSW_30nodes_120ranks.txt linear_with_manning 0.035 NSW > outfile.log
