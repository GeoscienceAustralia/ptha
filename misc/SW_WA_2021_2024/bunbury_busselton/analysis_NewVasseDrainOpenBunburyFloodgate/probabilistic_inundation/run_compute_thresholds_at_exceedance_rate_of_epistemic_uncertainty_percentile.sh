#!/bin/bash
#PBS -P w85
#PBS -q normalsr
#PBS -l walltime=24:00:00
#PBS -lmem=500GB
#PBS -lncpus=104
#PBS -ljobfs=20GB
#PBS -l wd
#PBS -l storage=scratch/w85+gdata/w85

source R_431_NCI_modules.sh

PERCENTILE=0.84 # epistemic uncertainty percentile
EXRATE=0.0004 # events per year 
MINDEPTH=0.0 # Lower limit to uniroot search (for efficiency)
MAXDEPTH=10.0 # Upper limit to uniroot search (for efficiency)
DEPTHTOL=0.001 # Tolerance for uniroot search
OUTPUTDIR=highres_domains_depth_at_epistemic_uncertainty_84pc_subsam_1

# In my study the high-resolution domains have a sequence of indices between these.
MIN_DOMAIN_INDEX=424
MAX_DOMAIN_INDEX=526

for domain_index in $(seq $MIN_DOMAIN_INDEX $MAX_DOMAIN_INDEX); do
    Rscript compute_threshold_at_exceedance_rate_of_epistemic_uncertainty_percentile.R depth $domain_index $PERCENTILE $EXRATE $MINDEPTH $MAXDEPTH $DEPTHTOL $OUTPUTDIR
done
