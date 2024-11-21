#!/bin/bash
#PBS -P w85
#PBS -q normal
#PBS -l walltime=2:00:00
#PBS -lmem=192GB
#PBS -lncpus=48
#PBS -l wd
#PBS -l storage=scratch/w85+gdata/w85

source ~/R_400_NCI_modules.sh
Rscript tar_multidomain_dirs.R

