#!/bin/bash
#PBS -P w85
#PBS -q normal
#PBS -l walltime=6:00:00
#PBS -lmem=190GB
#PBS -lncpus=48
#PBS -ljobfs=20GB
#PBS -l wd
#PBS -l storage=scratch/w85+gdata/w85+gdata/fj6

source R_431_NCI_modules.sh

PERCENTILE=0.84 # epistemic uncertainty percentile
EXRATE=0.0004 # events per year 
MINSEARCH=adaptive_minimum # Lower limit to uniroot search (for efficiency)
MAXSEARCH=adaptive_maximum # Upper limit to uniroot search (for efficiency)
UNIROOTTOL=0.001 # Tolerance for uniroot search

# In my study the high-resolution domains have a sequence of indices between these.
MIN_DOMAIN_INDEX=__MIN_DOMAIN_INDEX__
MAX_DOMAIN_INDEX=__MAX_DOMAIN_INDEX__
VARNAME=__VARNAME__

OUTPUTDIR="kalbarri2onslow_highres_domains_"$VARNAME"_percentile_"$PERCENTILE"_exrate_"$EXRATE"_hazard"

for domain_index in $(seq $MIN_DOMAIN_INDEX $MAX_DOMAIN_INDEX); do
    Rscript compute_threshold_at_exceedance_rate_of_epistemic_uncertainty_percentile.R $VARNAME $domain_index $PERCENTILE $EXRATE $MINSEARCH $MAXSEARCH $UNIROOTTOL $OUTPUTDIR
done
