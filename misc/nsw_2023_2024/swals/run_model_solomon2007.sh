#!/bin/bash
#PBS -P w85
#PBS -q normalsr
#PBS -l walltime=10:00:00
#PBS -lmem=1000GB
#PBS -lncpus=208
#PBS -l wd
#PBS -l storage=scratch/w85+gdata/w85

#source SWALS_ifort_modules_2023_llvm_B.sh
source SWALS_ifort_modules_2023_B.sh

#OMP_NUM_THREADS=13 mpiexec -np 16 --map-by ppr:4:socket:PE=13 ./model_sapphirerapids \
#    '../sources/like_historic/solomon2007/solomon2007-batch2_variable_area_uniform_slip_11244_count_1.tif' \
#    'multidomain_design_control_NNL4_defaultres.nml' \
#    run_solomon2007_11244 full 0.0 > outfile.log

#OMP_NUM_THREADS=13 mpiexec -np 16 --map-by ppr:4:socket:PE=13 ./model_sapphirerapids \
mpiexec -np 16 --map-by numa:SPAN --bind-to numa -x OMP_NUM_THREADS=13 -x OMP_PROC_BIND=true ./model_sapphirerapids \
    '../sources/like_historic/solomon2007/solomon2007-batch2_variable_area_uniform_slip_11244_count_1.tif' \
    'multidomain_design_control_NNL4_1arcminoffshore.nml' \
    run_solomon2007_11244_1arcminoffshore full 0.0 > outfile.log
