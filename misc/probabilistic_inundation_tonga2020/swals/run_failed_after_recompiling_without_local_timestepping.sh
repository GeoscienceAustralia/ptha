#!/bin/bash
#PBS -P w85
#PBS -q normal
#PBS -l walltime=10:00:00
#PBS -lmem=190GB
#PBS -lncpus=48
#PBS -l wd
#PBS -l storage=scratch/n74

source SWALS_ifort_modules.sh
ulimit -s unlimited

# This is the only run that blew-up. Try turning off local timestepping in the compilation and rerunning as single grids
OMP_NUM_THREADS=48 OMP_PROC_BIND=true mpiexec -np 1 --map-by ppr:1:node:PE=48 ./model ../sources/random/scenario_initial_conditions/kermadectonga2_row_0043831_Mw_95_HS.tif 0 0.8 ptha18_tonga_MSL0.8/ptha18_random_scenarios_kermadectonga2_row_0043831_Mw_95_HS 'full' '' linear_with_manning 0.035 tonga regional few_grids > outfile.log
#OMP_NUM_THREADS=24 OMP_PROC_BIND=true mpiexec -np 2 --map-by ppr:1:socket:PE=24 ./model ../sources/random/scenario_initial_conditions/kermadectonga2_row_0043831_Mw_95_HS.tif 0 0.8 ptha18_tonga_MSL0.8/ptha18_random_scenarios_kermadectonga2_row_0043831_Mw_95_HS 'full' '' linear_with_manning 0.035 tonga regional few_grids > outfile.log
