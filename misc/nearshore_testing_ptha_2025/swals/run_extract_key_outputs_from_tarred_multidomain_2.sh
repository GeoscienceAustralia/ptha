#!/bin/bash
#PBS -P w85
#PBS -q normal
#PBS -l walltime=6:00:00
#PBS -lmem=192GB
#PBS -lncpus=48
#PBS -l wd
#PBS -l storage=scratch/w85+gdata/w85

source R_431_NCI_modules.sh
Rscript extract_key_outputs_from_tarred_multidomain_2.R

