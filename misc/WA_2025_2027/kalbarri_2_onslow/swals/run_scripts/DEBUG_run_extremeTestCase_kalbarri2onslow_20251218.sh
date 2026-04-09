#!/bin/bash
#PBS -P w85
#PBS -q normalsr
#PBS -l walltime=2:00:00
#PBS -lmem=3000GB
#PBS -lncpus=624
#PBS -l wd
#PBS -l storage=scratch/w85+gdata/w85

source SWALS_ifx_modules_2025_llvm.sh

# Initial stage from Java 2006
stage_file=../sources/ptha18_test_sources/testing/VAUS_Mw92_sunda_arc/five_times_10231_for_extreme_testing.tif
# Default res and load balance. No tidal correction
multidomain_design_namelist=./multidomain_design_control/multidomain_kalbarri2onslow_20251218.nml
# Model name
model_name=extremeTestCase_kalbarri2onslow_20251218
# Ambient sea level of 0.0m AHD, approximation of expected sea level
ambient_sea_level=0

# Run the model
mpiexec -n 48 -map-by numa:SPAN -bind-to numa -x OMP_NUM_THREADS=13 -x OMP_PROC_BIND=true \
    ./model_debug \
    $stage_file \
    $multidomain_design_namelist \
    $model_name \
    full \
    $ambient_sea_level > ./log/$model_name.log$PBS_JOBID
