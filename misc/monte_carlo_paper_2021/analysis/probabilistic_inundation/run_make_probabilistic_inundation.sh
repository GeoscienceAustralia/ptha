#!/bin/bash
#PBS -P w85
#PBS -q normal
#PBS -l walltime=03:00:00
#PBS -lmem=190GB
#PBS -lncpus=48
#PBS -l wd
#PBS -l storage=gdata/n74+gdata/w85+gdata/fj6

source ~/R_400_NCI_modules.sh
source make_probabilistic_inundation.sh
