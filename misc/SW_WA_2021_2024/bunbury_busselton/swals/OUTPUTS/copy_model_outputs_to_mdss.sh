#!/bin/bash
#PBS -P w85
#PBS -q copyq
#PBS -l walltime=10:00:00
#PBS -lmem=4GB
#PBS -lncpus=1
#PBS -l wd
#PBS -l storage=scratch/w85+gdata/w85


# Copy the OUTPUTS folder for all scenarios to mdss
mdss put -r ptha18-BunburyBusseltonRevised-sealevel60cm tsunami/MODELS/inundation/WA_tsunami_inundation_DFES/bunbury_busselton/swals/OUTPUTS/
mdss put -r ptha18-BunburyBusseltonShutFloodgateRevised-sealevel60cm tsunami/MODELS/inundation/WA_tsunami_inundation_DFES/bunbury_busselton/swals/OUTPUTS/
