#!/bin/bash
#PBS -P w85
#PBS -q normal
#PBS -l walltime=48:00:00
#PBS -lmem=190GB
#PBS -lncpus=48
#PBS -ljobfs=50GB
#PBS -l wd
#PBS -l storage=scratch/w85+gdata/w85

source ../R_431_NCI_modules.sh

Rscript bzip2_some_files.R
