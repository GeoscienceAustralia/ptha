
#!/bin/bash
#PBS -P w85
#PBS -q normal
#PBS -l walltime=5:00:00
#PBS -lmem=1140GB
#PBS -lncpus=288
#PBS -l wd
#PBS -l storage=scratch/w85+gdata/w85

source SWALS_ifort_modules_2021.sh

OMP_NUM_THREADS=12 mpiexec -np 24 --map-by ppr:2:socket:PE=12 ./model '../sources/testing/VAUS_Mw92_sunda_arc/Mw92_scenarios_similar_size/Mw92_initial_condition_11144.tif' VAUS_Mw92_sunda_arc_Mw92_initial_condition_11144 full '../multidomain_design/domains_0.5_0.125/first_level_nesting.csv' '../multidomain_design/domains_0.5_0.125/second_level_nesting.csv' 'load_balance_files/load_balance_120821_9x_0.5_0.125_24MPI.txt' 0.0 > outfile.log

