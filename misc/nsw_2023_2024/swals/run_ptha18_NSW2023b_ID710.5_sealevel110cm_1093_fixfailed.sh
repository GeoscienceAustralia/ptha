#!/bin/bash
#PBS -P w85
#PBS -q normalsr
#PBS -l walltime=24:00:00
#PBS -lmem=1000GB
#PBS -lncpus=208
#PBS -l wd
#PBS -l storage=gdata/w85+scratch/w85

source SWALS_ifort_modules_2023_B.sh
# Load R as well (just for tarring directories)
module load R/4.3.1


mpiexec -np 16 --map-by numa:SPAN --bind-to numa -x OMP_NUM_THREADS=13 -x OMP_PROC_BIND=true ./model_sapphirerapids  \
    ../sources/hazard/scenarios_ID710.5/random_kermadectonga2/scenario_initial_conditions/kermadectonga2_row_0043461_Mw_94_HS.tif \
    multidomain_design_control_NNL4_1arcminoffshore.nml \
    ptha18-NSW2023b-ID710.5-sealevel110cm/random_kermadectonga2/ptha18_random_scenarios_kermadectonga2_row_0043461_Mw_94_HS \
    full \
    1.1 > outfile.log
# Tar the results 
Rscript tar_and_remove_matching_dir.R ./OUTPUTS/ptha18-NSW2023b-ID710.5-sealevel110cm/random_kermadectonga2/ptha18_random_scenarios_kermadectonga2_row_0043461_Mw_94_HS 

## mpiexec -np 16 --map-by numa:SPAN --bind-to numa -x OMP_NUM_THREADS=13 -x OMP_PROC_BIND=true ./model_sapphirerapids  \
##     ../sources/hazard/scenarios_ID710.5/random_outerrise_puysegur/scenario_initial_conditions/outerrise_puysegur_row_0003502_Mw_89_HS.tif \
##     multidomain_design_control_NNL4_1arcminoffshore.nml \
##     ptha18-NSW2023b-ID710.5-sealevel110cm/random_outerrise_puysegur/ptha18_random_scenarios_outerrise_puysegur_row_0003502_Mw_89_HS \
##     full \
##     1.1 > outfile.log
## # Tar the results 
## Rscript tar_and_remove_matching_dir.R ./OUTPUTS/ptha18-NSW2023b-ID710.5-sealevel110cm/random_outerrise_puysegur/ptha18_random_scenarios_outerrise_puysegur_row_0003502_Mw_89_HS 
## 
## mpiexec -np 16 --map-by numa:SPAN --bind-to numa -x OMP_NUM_THREADS=13 -x OMP_PROC_BIND=true ./model_sapphirerapids  \
##     ../sources/hazard/scenarios_ID710.5/random_southamerica/scenario_initial_conditions/southamerica_row_0147801_Mw_96_HS.tif \
##     multidomain_design_control_NNL4_1arcminoffshore.nml \
##     ptha18-NSW2023b-ID710.5-sealevel110cm/random_southamerica/ptha18_random_scenarios_southamerica_row_0147801_Mw_96_HS \
##     full \
##     1.1 > outfile.log
## # Tar the results 
## Rscript tar_and_remove_matching_dir.R ./OUTPUTS/ptha18-NSW2023b-ID710.5-sealevel110cm/random_southamerica/ptha18_random_scenarios_southamerica_row_0147801_Mw_96_HS 

