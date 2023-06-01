#!/bin/bash
#PBS -P w85
#PBS -q normal
#PBS -l walltime=02:00:00
#PBS -lmem=190GB
#PBS -lncpus=48
#PBS -ljobfs=20GB
#PBS -l wd
#PBS -l storage=scratch/w85+gdata/w85

source R_400_NCI_modules.sh

Rscript make_rasters.R
