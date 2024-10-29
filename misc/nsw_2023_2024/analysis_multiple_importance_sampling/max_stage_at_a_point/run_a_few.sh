#!/bin/bash
#PBS -P w85
#PBS -q normal
#PBS -l walltime=08:00:00
#PBS -lmem=96GB
#PBS -lncpus=24
#PBS -ljobfs=10GB
#PBS -l wd
#PBS -l storage=scratch/w85+gdata/w85+gdata/fj6

source ../../swals/R_431_NCI_modules.sh

for i in {1..4}; do
    Rscript extract_max_stage_at_a_point.R $i;
    done
