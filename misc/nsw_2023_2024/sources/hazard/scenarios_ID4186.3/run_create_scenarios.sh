#!/bin/bash
#PBS -P w85
#PBS -q normal
#PBS -l walltime=4:00:00
#PBS -l mem=190GB
#PBS -l ncpus=48
#PBS -l jobfs=20GB
#PBS -l wd
#PBS -l storage=scratch/w85+gdata/w85+gdata/fj6

source R_431_NCI_modules.sh

# Make the scenarios
Rscript create_scenarios.R

# Make the uplift/subsidence initial conditions
Rscript create_initial_conditions_for_scenarios.R

