#!/bin/bash
#PBS -P w85
#PBS -q normal
#PBS -l walltime=04:00:00
#PBS -lmem=16GB
#PBS -lncpus=4
#PBS -ljobfs=2GB
#PBS -l wd
#PBS -l storage=scratch/w85+gdata/w85+gdata/fj6

source ../../swals/R_431_NCI_modules.sh

for i in {1..4}; do
    Rscript extract_max_stage_at_a_point.R $i;
    done
