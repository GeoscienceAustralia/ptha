#!/bin/bash
#PBS -P w85
#PBS -q normal
#PBS -l walltime=06:00:00
#PBS -l mem=190GB
#PBS -l ncpus=48
#PBS -l jobfs=20GB
#PBS -l wd
#PBS -l storage=scratch/w85+gdata/w85

source ../../modules_R_431.sh

minrasts=_LOWER_
maxrasts=_UPPER_

# Deliberately skip domain 1 -- too large, kills memory
for i in $(seq $minrasts $maxrasts); do 
    Rscript compute_exceedance_rates_at_epistemic_uncertainty_percentile.R _VARIABLE_ _RANDOMSCENARIOS_ $i _PERCENTILE_ _EXCEEDANCETHRESHOLD_ _OUTPUTDIR_; 
done
