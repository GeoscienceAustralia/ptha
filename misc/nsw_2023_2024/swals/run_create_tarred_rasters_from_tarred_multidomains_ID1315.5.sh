#!/bin/bash
#PBS -P w85
#PBS -q normal
#PBS -l walltime=24:00:00
#PBS -lmem=190GB
#PBS -lncpus=48
#PBS -ljobfs=50GB
#PBS -l wd
#PBS -l storage=scratch/w85+gdata/w85

source R_431_NCI_modules.sh

# Split the work into this many separate jobs. Helps stop R accumulating memory and failing.
NRUNS=50

# Raster creation for scenario batch ID1315.5
for i in $(eval echo "{1..$NRUNS}"); do Rscript create_tarred_rasters_from_tarred_multidomains.R "OUTPUTS/ptha18-NSW2023-ID1315.5-sealevel110cm/random_*/ptha18_random_scenarios_*/RUN_*.tar" $i $NRUNS; done

