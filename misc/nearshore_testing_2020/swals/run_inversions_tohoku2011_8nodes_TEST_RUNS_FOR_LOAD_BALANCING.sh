#!/bin/bash
#PBS -P w85
#PBS -q normal
#PBS -l walltime=1:00:00
#PBS -lmem=1536GB
#PBS -lncpus=384
#PBS -l wd
#PBS -l storage=scratch/w85+gdata/w85

source SWALS_ifort_modules.sh
ulimit -s unlimited



#
# The load balancing will depend on the number of high-res areas (classified as NSW or Australia), and the offshore model type:
#     - linear_with_no_friction, linear_with_linear_friction can use the same load-balance info together
#     - linear_with_nonlinear_friction should use it's own load-balance info, because the offshore solve is more computationally demanding.
# Here we run models with a guess ('default')

# Australia-wide, linear -- short run that is used to create a 'real' load-balance file
OMP_NUM_THREADS=12 OMP_PROC_BIND=true mpiexec -np 32 --map-by ppr:2:socket:PE=12 ./model ../sources/Tohoku2011/RomanoEtAl2015/Tohoku2011_Romano_source_KAJIURA_SMOOTHED.tif 0.0 Tohoku2011_Romano2015 'test_load_balance' 'load_balance_files/load_balance_default_australia_8nodes_32mpi.txt' linear_with_linear_friction 0.0 australia > outfile.log

# NSW-only-high-res, linear  -- short run that is used to create a 'real' load-balance file
OMP_NUM_THREADS=12 OMP_PROC_BIND=true mpiexec -np 32 --map-by ppr:2:socket:PE=12 ./model ../sources/Tohoku2011/RomanoEtAl2015/Tohoku2011_Romano_source_KAJIURA_SMOOTHED.tif 0.0 Tohoku2011_Romano2015 'test_load_balance' 'load_balance_files/load_balance_default_NSW_8nodes_32mpi.txt' linear_with_linear_friction 0.0 NSW > outfile.log

# Australia-wide, linear+manning -- short run that is used to create a 'real' load-balance file
OMP_NUM_THREADS=12 OMP_PROC_BIND=true mpiexec -np 32 --map-by ppr:2:socket:PE=12 ./model ../sources/Tohoku2011/RomanoEtAl2015/Tohoku2011_Romano_source_KAJIURA_SMOOTHED.tif 0.0 Tohoku2011_Romano2015 'test_load_balance' 'load_balance_files/load_balance_default_australia_8nodes_32mpi.txt' linear_with_manning 0.035 australia > outfile.log

# NSW-wide, linear+manning -- short run that is used to create a 'real' load-balance file
OMP_NUM_THREADS=12 OMP_PROC_BIND=true mpiexec -np 32 --map-by ppr:2:socket:PE=12 ./model ../sources/Tohoku2011/RomanoEtAl2015/Tohoku2011_Romano_source_KAJIURA_SMOOTHED.tif 0.0 Tohoku2011_Romano2015 'test_load_balance' 'load_balance_files/load_balance_default_NSW_8nodes_32mpi.txt' linear_with_manning 0.035 NSW > outfile.log
