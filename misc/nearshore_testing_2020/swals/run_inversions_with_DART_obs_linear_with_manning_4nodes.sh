#!/bin/bash
#PBS -P w85
#PBS -q normal
#PBS -l walltime=18:00:00
#PBS -lmem=768GB
#PBS -lncpus=192
#PBS -l wd
#PBS -l storage=scratch/w85+gdata/w85

source SWALS_ifort_modules.sh
ulimit -s unlimited

# Choose number of openmp threads, such that 'mpi-processes'x'OMP_NUM_THREADS'='number of cpus'

#
# Here we run linear, global-only models for all source-inversions with DART data.
# (Chile 2010 (2 inversions); Tohoku 2011 (3 inversions); Chile 2015 (2 inversions)].
# The models have a manning-type friction term.
#

## Tohoku 2011 -- Romano et al 2014 
OMP_NUM_THREADS=8 OMP_PROC_BIND=true mpiexec -np 24 --map-by ppr:3:socket:PE=8 ./model ../sources/Tohoku2011/RomanoEtAl2015/Tohoku2011_Romano_source_KAJIURA_SMOOTHED.tif 0.0 Tohoku2011_Romano2015 full '' linear_with_manning 0.035 none > outfile.log
#
## Tohoku 2011 -- Satake et al 2013
OMP_NUM_THREADS=8 OMP_PROC_BIND=true mpiexec -np 24 --map-by ppr:3:socket:PE=8 ./model ../sources/Tohoku2011/SatakeEtAl2013/Satake_Tohoku11_sources_SUM_KAJIURA_SMOOTHED.tif  0.0 Tohoku2011_Satake2013 full '' linear_with_manning 0.035 none > outfile.log
#
## Tohoku 2011 -- Yamakazi et al 2018
OMP_NUM_THREADS=8 OMP_PROC_BIND=true mpiexec -np 24 --map-by ppr:3:socket:PE=8 ./model ../sources/Tohoku2011/YamakaziEtAl2018/yamakazi18_Tohoku11_sources_SUM_KAJIURA_SMOOTHED.tif 0.0 Tohoku2011_Yamakazi2018 full '' linear_with_manning 0.035 none > outfile.log
#
#
## Chile 2010 -- Fujii and Satake (2013)
OMP_NUM_THREADS=8 OMP_PROC_BIND=true mpiexec -np 24 --map-by ppr:3:socket:PE=8 ./model ../sources/Chile2010/FujiSatake2013/Fuji_chile2010_sources_SUM_KAJIURA_SMOOTHED.tif 0.0 Chile2010_Fujii2013 full '' linear_with_manning 0.035 none > outfile.log
#
## Chile 2010 -- Lorito et al (2011)
OMP_NUM_THREADS=8 OMP_PROC_BIND=true mpiexec -np 24 --map-by ppr:3:socket:PE=8 ./model ../sources/Chile2010/LoritoEtAl2011/Lorito_chile2010_sources_SUM_KAJIURA_SMOOTHED.tif 0.0 Chile2010_Lorito2011 full '' linear_with_manning 0.035 none > outfile.log
#
#
## Chile 2015 -- Williamson et al (2017)
OMP_NUM_THREADS=8 OMP_PROC_BIND=true mpiexec -np 24 --map-by ppr:3:socket:PE=8 ./model ../sources/Chile2015/WilliamsonEtAl2017/Williamson_chile_sources_SUM_KAJIURA_SMOOTHED.tif 0.0 Chile2015_Williamson2017 full '' linear_with_manning 0.035 none > outfile.log
#
## Chile 2015 -- Romano et al 2016
OMP_NUM_THREADS=8 OMP_PROC_BIND=true mpiexec -np 24 --map-by ppr:3:socket:PE=8 ./model ../sources/Chile2015/RomanoEtAl2016/Illapel_2015_Romano_KAJIURA_SMOOTHED.tif 0.0 Chile2015_Romano2016 full '' linear_with_manning 0.035 none > outfile.log

