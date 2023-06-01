#!/bin/bash
#PBS -P w85
#PBS -q normal
#PBS -l walltime=04:00:00
#PBS -lmem=190GB
#PBS -lncpus=48
#PBS -ljobfs=20GB
#PBS -l wd
#PBS -l storage=scratch/w85+gdata/w85

source R_400_NCI_modules.sh

minrasts=_LOWER_
maxrasts=_UPPER_

# Deliberately skip domain 1 -- too large, kills memory
for i in $(seq $minrasts $maxrasts); do 
    # Sunda calculation
    Rscript compute_exceedance_rates_at_epistemic_uncertainty_percentile.R depth ../../swals/OUTPUTS/_RUNDIR_/random_sunda2/ $i _PERCENTILE_ 0.001 _RUNDIR_-depth_exrate_0.001__PERCENTILE__sunda2; 
    # Outerrisesunda calculation
    Rscript compute_exceedance_rates_at_epistemic_uncertainty_percentile.R depth ../../swals/OUTPUTS/_RUNDIR_/random_outerrisesunda/ $i _PERCENTILE_ 0.001 _RUNDIR_-depth_exrate_0.001__PERCENTILE__outerrisesunda; 
done
