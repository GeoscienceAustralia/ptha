#!/bin/bash
#PBS -P w85
#PBS -q normal
#PBS -l walltime=4:00:00
#PBS -lmem=190GB
#PBS -lncpus=48
#PBS -ljobfs=20GB
#PBS -l wd
#PBS -l storage=scratch/w85+gdata/w85+gdata/fj6

source R_431_NCI_modules.sh

minrasts=_LOWER_
maxrasts=_UPPER_

# Deliberately skip domain 1 -- too large, kills memory
for i in $(seq $minrasts $maxrasts); do 
    Rscript compute_exceedance_rates_at_epistemic_uncertainty_percentile.R _VARIABLE_ _RANDOMSCENARIOS_ $i _PERCENTILE_ _EXCEEDANCETHRESHOLD_ _OUTPUTDIR_; 
done
