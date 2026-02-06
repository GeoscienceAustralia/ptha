#!/bin/bash
#PBS -P w85
#PBS -q copyq
#PBS -l walltime=10:00:00
#PBS -lmem=4GB
#PBS -lncpus=1
#PBS -l wd
#PBS -l storage=scratch/w85+gdata/w85


# Copy the OUTPUTS folder for all scenarios to mdss
#mdss mkdir tsunami/MODELS/inundation/australia_wide/clean_version/swals

mdss put -r OUTPUTS tsunami/MODELS/inundation/australia_wide/clean_version/swals/
mdss put -r OUTPUTS_SCRATCH tsunami/MODELS/inundation/australia_wide/clean_version/swals/
mdss put -r OUTPUTS_2022_new_events tsunami/MODELS/inundation/australia_wide/clean_version/swals/
mdss put -r OUTPUTS_2025_extend_WA tsunami/MODELS/inundation/australia_wide/clean_version/swals/
mdss put -r OUTPUTS_2025_NWWA tsunami/MODELS/inundation/australia_wide/clean_version/swals/
mdss put -r OUTPUTS_new_validation_events tsunami/MODELS/inundation/australia_wide/clean_version/swals/

