#!/bin/bash
#PBS -P n74
#PBS -q normal
#PBS -l walltime=12:00:00
#PBS -lmem=190GB
#PBS -lncpus=48
#PBS -l wd
#PBS -l storage=gdata/n74+gdata/w85

source ~/R_400_NCI_modules.sh
source make_depth_epistemic_uncertainty_rasters_16pc.sh
