#!/bin/bash
#PBS -P w85
#PBS -q normalsr
#PBS -l walltime=08:00:00
#PBS -lmem=500GB
#PBS -lncpus=104
#PBS -ljobfs=20GB
#PBS -l wd
#PBS -l storage=scratch/w85+gdata/w85

source R_431_NCI_modules.sh

Rscript compute_arrival_time_minima_and_scenario_average.R
