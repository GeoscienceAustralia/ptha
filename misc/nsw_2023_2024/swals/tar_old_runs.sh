#!/bin/bash
#PBS -P w85
#PBS -q normal
#PBS -l walltime=3:00:00
#PBS -lmem=96GB
#PBS -lncpus=24
#PBS -ljobfs=20GB
#PBS -l wd
#PBS -l storage=scratch/w85+gdata/w85

source R_431_NCI_modules.sh

for i in OUTPUTS/runs_before_Nov28_2023/run_*; do echo $i; Rscript tar_and_remove_matching_dir.R $i & done
wait
