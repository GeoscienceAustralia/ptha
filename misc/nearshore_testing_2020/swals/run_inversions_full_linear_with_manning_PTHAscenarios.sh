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

# Tohoku PTHA scenarios
OMP_NUM_THREADS=12 OMP_PROC_BIND=true mpiexec -np 32 --map-by ppr:2:socket:PE=12 ./model ../../ptha18_scenario_initial_conditions/kurilsjapan_46994_gauge_summary_stats_session_kurilsjapan_tohoku_2011_03_11_Mw9.1.Rdata.tif 0 Tohoku2011_ptha46994 full load_balance_files/load_balance_manning_offshore_NSW_8nodes_32ranks.txt linear_with_manning 0.035 NSW > outfile.log

OMP_NUM_THREADS=12 OMP_PROC_BIND=true mpiexec -np 32 --map-by ppr:2:socket:PE=12 ./model ../../ptha18_scenario_initial_conditions/kurilsjapan_47004_gauge_summary_stats_session_kurilsjapan_tohoku_2011_03_11_Mw9.1.Rdata.tif 0 Tohoku2011_ptha47004 full load_balance_files/load_balance_manning_offshore_NSW_8nodes_32ranks.txt linear_with_manning 0.035 NSW > outfile.log

OMP_NUM_THREADS=12 OMP_PROC_BIND=true mpiexec -np 32 --map-by ppr:2:socket:PE=12 ./model ../../ptha18_scenario_initial_conditions/kurilsjapan_47634_gauge_summary_stats_session_kurilsjapan_tohoku_2011_03_11_Mw9.1.Rdata.tif 0 Tohoku2011_ptha47634 full load_balance_files/load_balance_manning_offshore_NSW_8nodes_32ranks.txt linear_with_manning 0.035 NSW > outfile.log

# Chile 2010 PTHA scenarios
OMP_NUM_THREADS=12 OMP_PROC_BIND=true mpiexec -np 32 --map-by ppr:2:socket:PE=12 ./model ../../ptha18_scenario_initial_conditions/southamerica_128450_gauge_summary_stats_session_southamerica_chile_2010_02_27_Mw8.8.Rdata.tif 0 Chile2010_ptha128450 full load_balance_files/load_balance_manning_offshore_NSW_8nodes_32ranks.txt linear_with_manning 0.035 NSW > outfile.log

OMP_NUM_THREADS=12 OMP_PROC_BIND=true mpiexec -np 32 --map-by ppr:2:socket:PE=12 ./model ../../ptha18_scenario_initial_conditions/southamerica_128607_gauge_summary_stats_session_southamerica_chile_2010_02_27_Mw8.8.Rdata.tif 0 Chile2010_ptha128607 full load_balance_files/load_balance_manning_offshore_NSW_8nodes_32ranks.txt linear_with_manning 0.035 NSW > outfile.log

OMP_NUM_THREADS=12 OMP_PROC_BIND=true mpiexec -np 32 --map-by ppr:2:socket:PE=12 ./model ../../ptha18_scenario_initial_conditions/southamerica_128608_gauge_summary_stats_session_southamerica_chile_2010_02_27_Mw8.8.Rdata.tif 0 Chile2010_ptha128608 full load_balance_files/load_balance_manning_offshore_NSW_8nodes_32ranks.txt linear_with_manning 0.035 NSW > outfile.log

# Chile 2015 PTHA scenarios
OMP_NUM_THREADS=12 OMP_PROC_BIND=true mpiexec -np 32 --map-by ppr:2:socket:PE=12 ./model ../../ptha18_scenario_initial_conditions/southamerica_108269_gauge_summary_stats_session_southamerica_southamerica_2015_09_16_Mw8.3.Rdata.tif 0 Chile2015_ptha108269 full load_balance_files/load_balance_manning_offshore_NSW_8nodes_32ranks.txt linear_with_manning 0.035 NSW > outfile.log

OMP_NUM_THREADS=12 OMP_PROC_BIND=true mpiexec -np 32 --map-by ppr:2:socket:PE=12 ./model ../../ptha18_scenario_initial_conditions/southamerica_108275_gauge_summary_stats_session_southamerica_southamerica_2015_09_16_Mw8.3.Rdata.tif 0 Chile2015_ptha108275 full load_balance_files/load_balance_manning_offshore_NSW_8nodes_32ranks.txt linear_with_manning 0.035 NSW > outfile.log

OMP_NUM_THREADS=12 OMP_PROC_BIND=true mpiexec -np 32 --map-by ppr:2:socket:PE=12 ./model ../../ptha18_scenario_initial_conditions/southamerica_108314_gauge_summary_stats_session_southamerica_southamerica_2015_09_16_Mw8.3.Rdata.tif 0 Chile2015_ptha108314 full load_balance_files/load_balance_manning_offshore_NSW_8nodes_32ranks.txt linear_with_manning 0.035 NSW > outfile.log

# Chile 2014 PTHA scenarios
OMP_NUM_THREADS=12 OMP_PROC_BIND=true mpiexec -np 32 --map-by ppr:2:socket:PE=12 ./model ../../ptha18_scenario_initial_conditions/southamerica_87595_gauge_summary_stats_session_southamerica_southamerica_2014_04_01_Mw8.2.Rdata.tif 0 Chile2014_ptha87595 full load_balance_files/load_balance_manning_offshore_NSW_8nodes_32ranks.txt linear_with_manning 0.035 NSW > outfile.log

OMP_NUM_THREADS=12 OMP_PROC_BIND=true mpiexec -np 32 --map-by ppr:2:socket:PE=12 ./model ../../ptha18_scenario_initial_conditions/southamerica_87689_gauge_summary_stats_session_southamerica_southamerica_2014_04_01_Mw8.2.Rdata.tif 0 Chile2014_ptha87689 full load_balance_files/load_balance_manning_offshore_NSW_8nodes_32ranks.txt linear_with_manning 0.035 NSW > outfile.log

OMP_NUM_THREADS=12 OMP_PROC_BIND=true mpiexec -np 32 --map-by ppr:2:socket:PE=12 ./model ../../ptha18_scenario_initial_conditions/southamerica_87719_gauge_summary_stats_session_southamerica_southamerica_2014_04_01_Mw8.2.Rdata.tif 0 Chile2014_ptha87719 full load_balance_files/load_balance_manning_offshore_NSW_8nodes_32ranks.txt linear_with_manning 0.035 NSW > outfile.log

# Solomon 2007
OMP_NUM_THREADS=12 OMP_PROC_BIND=true mpiexec -np 32 --map-by ppr:2:socket:PE=12 ./model ../../ptha18_scenario_initial_conditions/solomon2_10634_gauge_summary_stats_session_solomon2_solomon_2007_04_01_Mw81.Rdata.tif 0 Solomon2007_ptha10634 full load_balance_files/load_balance_manning_offshore_NSW_8nodes_32ranks.txt linear_with_manning 0.035 NSW > outfile.log

OMP_NUM_THREADS=12 OMP_PROC_BIND=true mpiexec -np 32 --map-by ppr:2:socket:PE=12 ./model ../../ptha18_scenario_initial_conditions/solomon2_10637_gauge_summary_stats_session_solomon2_solomon_2007_04_01_Mw81.Rdata.tif 0 Solomon2007_ptha10637 full load_balance_files/load_balance_manning_offshore_NSW_8nodes_32ranks.txt linear_with_manning 0.035 NSW > outfile.log

OMP_NUM_THREADS=12 OMP_PROC_BIND=true mpiexec -np 32 --map-by ppr:2:socket:PE=12 ./model ../../ptha18_scenario_initial_conditions/solomon2_11216_gauge_summary_stats_session_solomon2_solomon_2007_04_01_Mw81.Rdata.tif 0 Solomon2007_ptha11216 full load_balance_files/load_balance_manning_offshore_NSW_8nodes_32ranks.txt linear_with_manning 0.035 NSW > outfile.log

# Puysegur 2009
OMP_NUM_THREADS=12 OMP_PROC_BIND=true mpiexec -np 32 --map-by ppr:2:socket:PE=12 ./model ../../ptha18_scenario_initial_conditions/puysegur2_1567_gauge_summary_stats_session_puysegur2_puysegur_2009_07_15_Mw7.8.Rdata.tif 0 Puysegur2009_ptha1567 full load_balance_files/load_balance_manning_offshore_NSW_8nodes_32ranks.txt linear_with_manning 0.035 NSW > outfile.log

OMP_NUM_THREADS=12 OMP_PROC_BIND=true mpiexec -np 32 --map-by ppr:2:socket:PE=12 ./model ../../ptha18_scenario_initial_conditions/puysegur2_1579_gauge_summary_stats_session_puysegur2_puysegur_2009_07_15_Mw7.8.Rdata.tif 0 Puysegur2009_ptha1579 full load_balance_files/load_balance_manning_offshore_NSW_8nodes_32ranks.txt linear_with_manning 0.035 NSW > outfile.log

OMP_NUM_THREADS=12 OMP_PROC_BIND=true mpiexec -np 32 --map-by ppr:2:socket:PE=12 ./model ../../ptha18_scenario_initial_conditions/puysegur2_1782_gauge_summary_stats_session_puysegur2_puysegur_2009_07_15_Mw7.8.Rdata.tif 0 Puysegur2009_ptha1782 full load_balance_files/load_balance_manning_offshore_NSW_8nodes_32ranks.txt linear_with_manning 0.035 NSW > outfile.log


