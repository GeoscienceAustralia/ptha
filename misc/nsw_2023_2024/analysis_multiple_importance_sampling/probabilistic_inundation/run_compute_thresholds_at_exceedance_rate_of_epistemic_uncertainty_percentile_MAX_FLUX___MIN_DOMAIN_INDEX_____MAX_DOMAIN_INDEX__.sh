#!/bin/bash
#PBS -P w85
#PBS -q normalsr
#PBS -l walltime=12:00:00
#PBS -lmem=500GB
#PBS -lncpus=104
#PBS -ljobfs=20GB
#PBS -l wd
#PBS -l storage=scratch/w85+gdata/w85+gdata/fj6

source R_431_NCI_modules.sh

PERCENTILE=0.84 # epistemic uncertainty percentile
EXRATE=0.0004 # events per year 
MINSTAGE=adaptive_minimum # Lower limit to uniroot search (for efficiency)
MAXSTAGE=adaptive_maximum # Upper limit to uniroot search (for efficiency)
STAGETOL=0.01 # Tolerance for uniroot search
OUTPUTDIR=nsw_full_coast_MIS_highres_domains_max_flux_at_epistemic_uncertainty_84pc

# In my study the high-resolution domains have a sequence of indices between these.
MIN_DOMAIN_INDEX=__MIN_DOMAIN_INDEX__
MAX_DOMAIN_INDEX=__MAX_DOMAIN_INDEX__

for domain_index in $(seq $MIN_DOMAIN_INDEX $MAX_DOMAIN_INDEX); do
    Rscript compute_threshold_at_exceedance_rate_of_epistemic_uncertainty_percentile.R max_flux $domain_index $PERCENTILE $EXRATE $MINSTAGE $MAXSTAGE $STAGETOL $OUTPUTDIR
done
