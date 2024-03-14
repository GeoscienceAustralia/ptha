#!/bin/bash
#PBS -P w85
#PBS -q copyq
#PBS -l walltime=10:00:00
#PBS -lmem=4GB
#PBS -lncpus=1
#PBS -l wd
#PBS -l storage=scratch/w85+gdata/w85


# Copy the OUTPUTS folder for all scenarios to mdss
mdss mkdir tsunami/MODELS/inundation/NSW_tsunami_inundation_with_NSWSES/stage1/full_coast/swals/OUTPUTS/
mdss put -r ptha18-NSW2023b-ID710.5-sealevel110cm tsunami/MODELS/inundation/NSW_tsunami_inundation_with_NSWSES/stage1/full_coast/swals/OUTPUTS/
