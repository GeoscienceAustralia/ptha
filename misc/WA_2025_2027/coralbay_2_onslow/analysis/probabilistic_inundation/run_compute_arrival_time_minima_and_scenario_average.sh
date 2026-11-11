#!/bin/bash
#PBS -P w85
#PBS -q normal
#PBS -l walltime=08:00:00
#PBS -l mem=190GB
#PBS -l ncpus=48
#PBS -l jobfs=20GB
#PBS -l wd
#PBS -l storage=scratch/w85+gdata/w85

source ../../modules_R_431.sh

Rscript compute_arrival_time_minima_and_scenario_average.R
