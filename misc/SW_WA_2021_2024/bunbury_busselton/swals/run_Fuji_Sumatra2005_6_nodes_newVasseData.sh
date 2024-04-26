#!/bin/bash
#PBS -P w85
#PBS -q normal
#PBS -l walltime=04:00:00
#PBS -lmem=1140GB
#PBS -lncpus=288
#PBS -l wd
#PBS -l storage=scratch/w85+gdata/w85

source SWALS_ifort_modules_2023.sh

# Full model
OMP_NUM_THREADS=12 mpiexec -np 24 --map-by ppr:2:socket:PE=12 ./model \
    '../sources/like_historic/sumatra2005/Fuji_nias2005_unit_sources_SUM_KAJIURA_SMOOTHED.tif' \
    Fuji_sumatra2005_newVasseData_24hrs_domain301122 \
    full \
    "../multidomain_design/domains_301122_0.5_0.166666666666667_0.0333333333333333_0.00555555555555556/first_level_nesting_edited.csv" \
    "../multidomain_design/domains_301122_0.5_0.166666666666667_0.0333333333333333_0.00555555555555556/second_level_nesting_edited.csv" \
    "../multidomain_design/domains_301122_0.5_0.166666666666667_0.0333333333333333_0.00555555555555556/third_level_nesting_edited.csv" \
    "../multidomain_design/domains_301122_0.5_0.166666666666667_0.0333333333333333_0.00555555555555556/fourth_level_nesting_edited.csv" \
    "./load_balance_files/load_balance_301122_24MPIb.txt" \
    0.0 > outfile.log
