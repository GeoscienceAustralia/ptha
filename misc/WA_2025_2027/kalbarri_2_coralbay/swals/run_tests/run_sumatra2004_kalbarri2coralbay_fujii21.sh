#!/bin/bash
#PBS -P w85
#PBS -q normalsr
#PBS -l walltime=6:00:00
#PBS -lmem=1000GB
#PBS -lncpus=208
#PBS -l wd
#PBS -l storage=scratch/w85+gdata/w85

source SWALS_ifx_modules_2025_llvm.sh

stage_file=../sources/like_historic/sumatra2004/Fujii2021/Fujii21_time_varying_forcing_realistic.csv
## No tidal adjustment
multidomain_design_namelist=./multidomain_design_control/multidomain_kalbarri2coralbay_B_historicaltsunamis.nml

# Model name
model_name=sumatra2004_Fujii21_timevarying
# Ambient sea level of 0.0m AHD, approximation of expected sea level
ambient_sea_level=0

# Run the model
mpiexec -n 16 -map-by numa:SPAN -bind-to numa -x OMP_NUM_THREADS=13 -x OMP_PROC_BIND=true \
    ./model_build \
    $stage_file \
    $multidomain_design_namelist \
    $model_name \
    full \
    $ambient_sea_level > outfile.log
