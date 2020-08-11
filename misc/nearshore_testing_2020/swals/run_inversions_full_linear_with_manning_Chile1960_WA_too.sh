#!/bin/bash
#PBS -P w85
#PBS -q normal
#PBS -l walltime=12:00:00
#PBS -lmem=1536GB
#PBS -lncpus=384
#PBS -l wd
#PBS -l storage=scratch/w85+gdata/w85

source SWALS_ifort_modules.sh
ulimit -s unlimited

OMP_NUM_THREADS=12 OMP_PROC_BIND=true mpiexec -np 32 --map-by ppr:2:socket:PE=12 ./model ../sources/Chile1960/FujiSatake2013/Fuji_chile1960_sources_SUM_KAJIURA_SMOOTHED.tif 0 Chile1960_FujiSatake2013 full load_balance_files/load_balance_manning_offshore_australia_8nodes_32ranks.txt linear_with_manning 0.035 australia > outfile.log
OMP_NUM_THREADS=12 OMP_PROC_BIND=true mpiexec -np 32 --map-by ppr:2:socket:PE=12 ./model ../sources/Chile1960/HoEtAl2019/Ho_Chile1960_initial_displacement.tif 0 Chile1960_HoEtAl2019 full load_balance_files/load_balance_manning_offshore_australia_8nodes_32ranks.txt linear_with_manning 0.035 australia > outfile.log
