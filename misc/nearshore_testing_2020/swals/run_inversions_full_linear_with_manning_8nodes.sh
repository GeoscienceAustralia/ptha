#!/bin/bash
#PBS -P w85
#PBS -q normal
#PBS -l walltime=48:00:00
#PBS -lmem=1536GB
#PBS -lncpus=384
#PBS -l wd
#PBS -l storage=scratch/w85+gdata/w85

source SWALS_ifort_modules.sh
ulimit -s unlimited

OMP_NUM_THREADS=12 OMP_PROC_BIND=true mpiexec -np 32 --map-by ppr:2:socket:PE=12 ./model ../sources/Chile1960/FujiSatake2013/Fuji_chile1960_sources_SUM_KAJIURA_SMOOTHED.tif 0 Chile1960_FujiSatake2013 full load_balance_files/load_balance_manning_offshore_NSW_8nodes_32ranks.txt linear_with_manning 0.035 NSW > outfile.log
OMP_NUM_THREADS=12 OMP_PROC_BIND=true mpiexec -np 32 --map-by ppr:2:socket:PE=12 ./model ../sources/Chile1960/HoEtAl2019/Ho_Chile1960_initial_displacement.tif 0 Chile1960_HoEtAl2019 full load_balance_files/load_balance_manning_offshore_NSW_8nodes_32ranks.txt linear_with_manning 0.035 NSW > outfile.log
OMP_NUM_THREADS=12 OMP_PROC_BIND=true mpiexec -np 32 --map-by ppr:2:socket:PE=12 ./model ../sources/Sumatra2004/FujiSatake2007/Fuji_andaman2004_unit_sources_SUM_KAJIURA_SMOOTHED.tif 0 Sumatra2004_FujiSatake2007 full load_balance_files/load_balance_manning_offshore_australia_8nodes_32ranks.txt linear_with_manning 0.035 australia > outfile.log
OMP_NUM_THREADS=12 OMP_PROC_BIND=true mpiexec -np 32 --map-by ppr:2:socket:PE=12 ./model ../sources/Sumatra2004/LoritoEtAl2010/Lorito_sumatra2004_sources_SUM_KAJIURA_SMOOTHED.tif 0 Sumatra2004_LoritoEtAl2010 full load_balance_files/load_balance_manning_offshore_australia_8nodes_32ranks.txt linear_with_manning 0.035 australia > outfile.log
OMP_NUM_THREADS=12 OMP_PROC_BIND=true mpiexec -np 32 --map-by ppr:2:socket:PE=12 ./model ../sources/Sumatra2004/PiatanesiLorito2007/Piatanesi_sumatra2004_sources_SUM_KAJIURA_SMOOTHED.tif 0 Sumatra2004_PiatanesiLorito2007 full load_balance_files/load_balance_manning_offshore_australia_8nodes_32ranks.txt linear_with_manning 0.035 australia > outfile.log
OMP_NUM_THREADS=12 OMP_PROC_BIND=true mpiexec -np 32 --map-by ppr:2:socket:PE=12 ./model ../sources/Chile2010/FujiSatake2013/Fuji_chile2010_sources_SUM_KAJIURA_SMOOTHED.tif 0 Chile2010_FujiSatake2013 full load_balance_files/load_balance_manning_offshore_NSW_8nodes_32ranks.txt linear_with_manning 0.035 NSW > outfile.log
OMP_NUM_THREADS=12 OMP_PROC_BIND=true mpiexec -np 32 --map-by ppr:2:socket:PE=12 ./model ../sources/Chile2010/LoritoEtAl2011/Lorito_chile2010_sources_SUM_KAJIURA_SMOOTHED.tif 0 Chile2010_LoritoEtAl2011 full load_balance_files/load_balance_manning_offshore_NSW_8nodes_32ranks.txt linear_with_manning 0.035 NSW > outfile.log
OMP_NUM_THREADS=12 OMP_PROC_BIND=true mpiexec -np 32 --map-by ppr:2:socket:PE=12 ./model ../sources/Tohoku2011/SatakeEtAl2013/Satake_Tohoku11_sources_SUM_KAJIURA_SMOOTHED.tif 0 Tohoku2011_SatakeEtAl2013 full load_balance_files/load_balance_manning_offshore_NSW_8nodes_32ranks.txt linear_with_manning 0.035 NSW > outfile.log
OMP_NUM_THREADS=12 OMP_PROC_BIND=true mpiexec -np 32 --map-by ppr:2:socket:PE=12 ./model ../sources/Tohoku2011/YamakaziEtAl2018/yamakazi18_Tohoku11_sources_SUM_KAJIURA_SMOOTHED.tif 0 Tohoku2011_YamakaziEtAl2018 full load_balance_files/load_balance_manning_offshore_NSW_8nodes_32ranks.txt linear_with_manning 0.035 NSW > outfile.log
OMP_NUM_THREADS=12 OMP_PROC_BIND=true mpiexec -np 32 --map-by ppr:2:socket:PE=12 ./model ../sources/Tohoku2011/RomanoEtAl2015/Tohoku2011_Romano_source_KAJIURA_SMOOTHED.tif 0 Tohoku2011_RomanoEtAl2015 full load_balance_files/load_balance_manning_offshore_NSW_8nodes_32ranks.txt linear_with_manning 0.035 NSW > outfile.log
OMP_NUM_THREADS=12 OMP_PROC_BIND=true mpiexec -np 32 --map-by ppr:2:socket:PE=12 ./model ../sources/Chile2015/WilliamsonEtAl2017/Williamson_chile_sources_SUM_KAJIURA_SMOOTHED.tif 0 Chile2015_WilliamsonEtAl2017 full load_balance_files/load_balance_manning_offshore_NSW_8nodes_32ranks.txt linear_with_manning 0.035 NSW > outfile.log
OMP_NUM_THREADS=12 OMP_PROC_BIND=true mpiexec -np 32 --map-by ppr:2:socket:PE=12 ./model ../sources/Chile2015/RomanoEtAl2016/Illapel_2015_Romano_KAJIURA_SMOOTHED.tif 0 Chile2015_RomanoEtAl2016 full load_balance_files/load_balance_manning_offshore_NSW_8nodes_32ranks.txt linear_with_manning 0.035 NSW > outfile.log

