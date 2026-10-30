#!/bin/bash
#PBS -P w85
#PBS -q normal
#PBS -l walltime=15:00:00
#PBS -l mem=190GB
#PBS -l ncpus=48
#PBS -l jobfs=20GB
#PBS -l wd
#PBS -l storage=scratch/w85+gdata/w85+gdata/fj6

source ../../modules_R_431.sh

PERCENTILE=0.84 # epistemic uncertainty percentile
EXRATE=0.0004 # events per year 
MINSEARCH=adaptive_minimum # Lower limit to uniroot search (for efficiency)
MAXSEARCH=adaptive_maximum # Upper limit to uniroot search (for efficiency)
UNIROOTTOL=0.001 # Tolerance for uniroot search

# In my study the high-resolution domains have a sequence of indices between these.
MIN_DOMAIN_INDEX=231
MAX_DOMAIN_INDEX=239
VARNAME=max_flux

OUTPUTDIR="coralbay2onslow_highres_domains_"$VARNAME"_percentile_"$PERCENTILE"_exrate_"$EXRATE"_hazard"

#FIXME: Put the problematic domain indices here
#for domain_index in 230 239 247 248; do
#    Rscript compute_threshold_at_exceedance_rate_of_epistemic_uncertainty_percentile.R $VARNAME $domain_index $PERCENTILE $EXRATE $MINSEARCH $MAXSEARCH $UNIROOTTOL $OUTPUTDIR
#done
