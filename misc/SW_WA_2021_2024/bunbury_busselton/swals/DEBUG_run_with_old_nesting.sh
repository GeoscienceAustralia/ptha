#!/bin/bash
#PBS -P w85
#PBS -q normal
#PBS -l walltime=04:00:00
#PBS -lmem=1140GB
#PBS -lncpus=288
#PBS -l wd
#PBS -l storage=gdata/w85+scratch/w85

source SWALS_ifort_modules_2023.sh
# Load R as well (just for tarring directories)
module load R/4.2.1 

#
# This model went unstable with the "new" SWALS nesting (the only scenario that did)
# I re-ran with -DOLD_PROCESS_DATA_TO_SEND_B4FEB22
#

OMP_NUM_THREADS=12 OMP_PROC_BIND=true mpiexec -np 24 --map-by ppr:2:socket:PE=12 ./model_OLD_NESTING  \
    ../sources/hazard/random_sunda2/scenario_initial_conditions/sunda2_row_0107982_Mw_94_HS.tif \
    ptha18-BunburyBusseltonRevised-sealevel60cm/random_sunda2/ptha18_random_scenarios_sunda2_row_0107982_Mw_94_HS \
    full \
    ../multidomain_design/domains_301122_0.5_0.166666666666667_0.0333333333333333_0.00555555555555556/first_level_nesting_edited.csv \
    ../multidomain_design/domains_301122_0.5_0.166666666666667_0.0333333333333333_0.00555555555555556/second_level_nesting_edited.csv \
    ../multidomain_design/domains_301122_0.5_0.166666666666667_0.0333333333333333_0.00555555555555556/third_level_nesting_edited.csv \
    ../multidomain_design/domains_301122_0.5_0.166666666666667_0.0333333333333333_0.00555555555555556/fourth_level_nesting_edited.csv \
    load_balance_files/load_balance_301122_24MPIb.txt \
    0.6 > outfile.log
# Tar the results 
Rscript tar_and_remove_matching_dir.R ./OUTPUTS/ptha18-BunburyBusseltonRevised-sealevel60cm/random_sunda2/ptha18_random_scenarios_sunda2_row_0107982_Mw_94_HS

