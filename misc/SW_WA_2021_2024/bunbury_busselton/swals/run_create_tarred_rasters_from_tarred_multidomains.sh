#!/bin/bash
#PBS -P w85
#PBS -q normal
#PBS -l walltime=10:00:00
#PBS -lmem=190GB
#PBS -lncpus=48
#PBS -ljobfs=20GB
#PBS -l wd
#PBS -l storage=scratch/w85+gdata/w85

source R_421_NCI_modules.sh

for i in {1..50}; do Rscript create_tarred_rasters_from_tarred_multidomains.R "OUTPUTS/ptha18-BunburyBusseltonRevised-sealevel60cm/random_*/ptha18_random_scenarios_*/RUN_*.tar" $i 50; done
