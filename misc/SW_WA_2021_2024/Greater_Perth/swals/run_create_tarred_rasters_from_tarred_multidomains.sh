#!/bin/bash
#PBS -P w85
#PBS -q normal
#PBS -l walltime=10:00:00
#PBS -lmem=190GB
#PBS -lncpus=48
#PBS -ljobfs=20GB
#PBS -l wd
#PBS -l storage=scratch/w85+gdata/w85

source R_400_NCI_modules.sh

#for i in {1..50}; do Rscript create_tarred_rasters_from_tarred_multidomains.R "OUTPUTS/ptha18-GreaterPerth-sealevel60cm-lowres/random_*/ptha18_random_scenarios_*/RUN_*.tar" $i 50; done

for i in {1..50}; do Rscript create_tarred_rasters_from_tarred_multidomains.R "OUTPUTS/ptha18-GreaterPerth-sealevel60cm-reviseddomain-highres/random_*/ptha18_random_scenarios_*/RUN_*.tar" $i 50; done

#for i in {1..1}; do Rscript create_tarred_rasters_from_tarred_multidomains.R "OUTPUTS/ptha18-GreaterPerth-sealevel60cm-*/random_*/ptha18_random_scenarios_*/RUN_*.tar" $i 1; done

#Rscript create_tarred_rasters_from_tarred_multidomains.R "OUTPUTS/ptha18-GreaterPerth-sealevel60cm/random_sunda2/ptha18_random_scenarios_sunda2_row_0110918_Mw_96_HS-full-ambient_sea_level_0.6/RUN*.tar" 1 1

#Rscript create_tarred_rasters_from_tarred_multidomains.R "OUTPUTS/ptha18-GreaterPerth-sealevel60cm/random_sunda2/ptha18_random_scenarios_sunda2_row_0108445_Mw_94_HS-full-ambient_sea_level_0.6/RUN*.tar" 1 1
