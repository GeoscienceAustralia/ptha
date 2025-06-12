#!/bin/bash
#PBS -P w85
#PBS -q normal
#PBS -l walltime=48:00:00
#PBS -lmem=1536GB
#PBS -lncpus=384
#PBS -l wd
#PBS -l storage=scratch/w85+gdata/w85

#
# This script runs the source inversions as used for the paper. In practice
# they were run with the same commands but in separate scripts -- this is
# provided to make clearer what was done.
#


source SWALS_ifort_modules_2023.sh

# Chile 1960
OMP_NUM_THREADS=12 OMP_PROC_BIND=true mpiexec -np 32 --map-by ppr:2:socket:PE=12 ./model ../sources/Chile1960/HoEtAl2019/Ho_Chile1960_initial_displacement.tif 0 Chile1960_HoEtAl2019 full load_balance_files/load_balance_manning_offshore_NSW_8nodes_32ranks.txt linear_with_manning 0.035 NSW > outfile.log

# Sumatra 2004 with a finite propagation speed
mpiexec -np 32 --map-by ppr:2:socket:PE=12 -x OMP_NUM_THREADS=12 -x OMP_PROC_BIND=true  ./model ../sources/sources_with_rise_time/Sumatra2004/FujiSatake2007/FujiSatake2007_time_varying_forcing_realistic.csv 0.0 Sumatra2004_FujiiandSatake2007_time_varying full load_balance_files/load_balance_australiaWA_8nodes_32mpi.txt linear_with_manning 0.035 australiaWA > outfile.log

# Sumatra 2005
mpiexec -np 32 --map-by ppr:2:socket:PE=12 -x OMP_NUM_THREADS=12 -x OMP_PROC_BIND=true  ./model ../sources/sumatra2005/Fuji_nias2005_unit_sources_SUM_KAJIURA_SMOOTHED.tif 0 Sumatra2005_Fujiietal2020 full load_balance_files/load_balance_wa_8nodes_32MPI.txt linear_with_manning 0.035 WA > outfile.log

# Java2006 with a finite propagation speed
mpiexec -np 32 --map-by ppr:2:socket:PE=12 -x OMP_NUM_THREADS=12 -x OMP_PROC_BIND=true  ./model ../sources/sources_with_rise_time/Java2006/FujiiSatake06/FujiSatake2006_Java2006_time_varying_forcing_realistic.csv 0.0 Java2006_FujiiandSatake2006_time_varying full load_balance_files/load_balance_nwwa_8nodes_32mpi.txt linear_with_manning 0.035 NWWA > outfile.log

# Solomon 2007 
mpiexec -np 32 --map-by ppr:2:socket:PE=12 -x OMP_NUM_THREADS=12 -x OMP_PROC_BIND=true  ./model ../sources/Solomon2007/Wei_2015/Wei_S2_Solomon2007_source_SUM_KAJIURA_SMOOTHED.tif 0 Solomon2007_Weietal2015 full load_balance_files/load_balance_manning_offshore_NSW_8nodes_32ranks.txt linear_with_manning 0.035 NSW > outfile.log

# Sumatra 2007
mpiexec -np 32 --map-by ppr:2:socket:PE=12 -x OMP_NUM_THREADS=12 -x OMP_PROC_BIND=true  ./model ../sources/Sumatra2007/FujiiSatake2007/FujiiSatake08_2007_Sumatra_SUM_KAJIURA_SMOOTHED.tif 0 Sumatra2007_Fujiietal2008 full load_balance_files/load_balance_nwwa_8nodes_32mpi.txt linear_with_manning 0.035 NWWA > outfile.log

# Puysegur 2009 
mpiexec -np 32 --map-by ppr:2:socket:PE=12 -x OMP_NUM_THREADS=12 -x OMP_PROC_BIND=true  ./model ../sources/Puysegur2009/Bevan2010/Bevan_Puysegur2009_source_SUM_KAJIURA_SMOOTHED.tif 0 Puysegur2009_Bevanetal2010 full load_balance_files/load_balance_manning_offshore_NSW_8nodes_32ranks.txt linear_with_manning 0.035 NSW > outfile.log

# Chile 2010
OMP_NUM_THREADS=12 OMP_PROC_BIND=true mpiexec -np 32 --map-by ppr:2:socket:PE=12 ./model ../sources/Chile2010/LoritoEtAl2011/Lorito_chile2010_sources_SUM_KAJIURA_SMOOTHED.tif 0 Chile2010_LoritoEtAl2011 full load_balance_files/load_balance_manning_offshore_NSW_8nodes_32ranks.txt linear_with_manning 0.035 NSW > outfile.log

# Tohoku 2011
OMP_NUM_THREADS=12 OMP_PROC_BIND=true mpiexec -np 32 --map-by ppr:2:socket:PE=12 ./model ../sources/Tohoku2011/YamazakiEtAl2018Fixed/yamazaki18_Tohoku11_sources_SUM_KAJIURA_SMOOTHED.tif 0 Tohoku2011_YamazakiEtAl2018Fixed full load_balance_files/load_balance_manning_offshore_NSW_8nodes_32ranks.txt linear_with_manning 0.035 NSW > outfile.log

# Chile 2014 
mpiexec -np 32 --map-by ppr:2:socket:PE=12 -x OMP_NUM_THREADS=12 -x OMP_PROC_BIND=true  ./model ../sources/Chile2014/An_2014/An_Chile2014_source_SUM_KAJIURA_SMOOTHED.tif 0 Chile2014_AnEtAl2014 full load_balance_files/load_balance_manning_offshore_NSW_8nodes_32ranks.txt linear_with_manning 0.035 NSW > outfile.log

# Chile 2015
OMP_NUM_THREADS=12 OMP_PROC_BIND=true mpiexec -np 32 --map-by ppr:2:socket:PE=12 ./model ../sources/Chile2015/WilliamsonEtAl2017/Williamson_chile_sources_SUM_KAJIURA_SMOOTHED.tif 0 Chile2015_WilliamsonEtAl2017 full load_balance_files/load_balance_manning_offshore_NSW_8nodes_32ranks.txt linear_with_manning 0.035 NSW > outfile.log

# Sandwich 2021
mpiexec -np 32 --map-by ppr:2:socket:PE=12 -x OMP_NUM_THREADS=12 -x OMP_PROC_BIND=true  ./model ../sources/Sandwich2021/Roger2024_GCMT_source/RogerGCMT_2024_Sandwich2021_SUM_KAJIURA_SMOOTHED.tif 0 Sandwich2021_Rogeretal2024 full load_balance_files/load_balance_manning_offshore_perth_8nodes_32ranks_2022.txt linear_with_manning 0.035 perth > outfile.log

# New Hebrides 2021
mpiexec -np 32 --map-by ppr:2:socket:PE=12 -x OMP_NUM_THREADS=12 -x OMP_PROC_BIND=true  ./model ../sources/NewHebrides2021/GusmanEtAl/Gusman_NewHebrides_sources_SUM_KAJIURA_SMOOTHED.tif 0 NewHebrides2021_GusmanEtAl2022 full load_balance_files/load_balance_manning_offshore_NSW_8nodes_32ranks.txt linear_with_manning 0.035 NSW > outfile.log

# Kermadec 2021
OMP_NUM_THREADS=12 OMP_PROC_BIND=true mpiexec -np 32 --map-by ppr:2:socket:PE=12 ./model ../sources/KermadecTonga2021/Romano2021/Kermadec2021_Romano.tif 0 Kermadec2021_Romano_source full load_balance_files/load_balance_manning_offshore_NSW_8nodes_32ranks_2022.txt linear_with_manning 0.035 NSW > outfile.log

