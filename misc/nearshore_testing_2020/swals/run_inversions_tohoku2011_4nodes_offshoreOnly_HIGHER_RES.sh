#!/bin/bash
#PBS -P w85
#PBS -q normal
#PBS -l walltime=10:00:00
#PBS -lmem=768GB
#PBS -lncpus=192
#PBS -l wd
#PBS -l storage=scratch/w85+gdata/w85

source SWALS_ifort_modules.sh
ulimit -s unlimited

# Choose number of openmp threads, such that 'mpi-processes'x'OMP_NUM_THREADS'='number of cpus'

#
# The run here uses mesh_refine=8 (i.e. 2x mesh refinement), and is related to the regular runs commented out below which
# used mesh_refine=4. The code was recompiled to support the higher mesh-refinement.
# The point is to confirm that the late-time waves are 'basically convergent'
# -- i.e. are not attributable to the use of an overly coarse grid or large timestep.
#
OMP_NUM_THREADS=4 OMP_PROC_BIND=true mpiexec -np 48 --map-by ppr:6:socket:PE=4 ./model ../sources/Tohoku2011/RomanoEtAl2015/Tohoku2011_Romano_source_KAJIURA_SMOOTHED.tif 0.0 Tohoku2011_Romano2015_HIGHER_RES full '' linear_with_no_friction 0.0 none > outfile.log


##
## These 3 runs use the standard setup -- they were used to check the importance of rise-time, using a relatively rough source-model.
##
#
#OMP_NUM_THREADS=4 OMP_PROC_BIND=true mpiexec -np 24 --map-by ppr:6:socket:PE=4 ./model ../sources/Tohoku2011/RomanoEtAl2015/Tohoku2011_Romano_source_KAJIURA_SMOOTHED.tif 0.0 Tohoku2011_Romano2015 full '' linear_with_no_friction 0.0 none > outfile.log
#
#OMP_NUM_THREADS=4 OMP_PROC_BIND=true mpiexec -np 24 --map-by ppr:6:socket:PE=4 ./model ../sources/Tohoku2011/RomanoEtAl2015/Tohoku2011_Romano_source_KAJIURA_SMOOTHED.tif 60.0 Tohoku2011_Romano2015 full '' linear_with_no_friction 0.0 none > outfile.log
#
#OMP_NUM_THREADS=4 OMP_PROC_BIND=true mpiexec -np 24 --map-by ppr:6:socket:PE=4 ./model ../sources/Tohoku2011/RomanoEtAl2015/Tohoku2011_Romano_source_KAJIURA_SMOOTHED.tif 30.0 Tohoku2011_Romano2015 full '' linear_with_no_friction 0.0 none > outfile.log

