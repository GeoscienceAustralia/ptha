#!/bin/bash
#PBS -P w85
#PBS -q normal
#PBS -l walltime=03:00:00
#PBS -lmem=190GB
#PBS -lncpus=48
#PBS -ljobfs=20GB
#PBS -l wd
#PBS -l storage=scratch/w85+gdata/w85

source R_400_NCI_modules.sh

Rscript create_plots_from_tarred_multidomain_dirs.R "OUTPUTS/ptha18-GreaterPerth-sealevel60cm-revised/random*/ptha*/RUN*.tar"
