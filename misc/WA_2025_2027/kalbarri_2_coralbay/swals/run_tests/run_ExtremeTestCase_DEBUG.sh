#!/bin/bash
#PBS -P w85
#PBS -q normalsr
#PBS -l walltime=12:00:00
#PBS -lmem=1000GB
#PBS -lncpus=208
#PBS -l wd
#PBS -l storage=scratch/w85+gdata/w85

source SWALS_ifx_modules_2025_llvm.sh

# Initial stage 
stage_file=../sources/test_sources/testing/VAUS_Mw92_sunda_arc/five_times_10231_for_extreme_testing.tif

## With tidal adjustment
multidomain_design_namelist=./multidomain_design_control/multidomain_kalbarri2coralbay_B_hazard.nml

## No tidal adjustment
#multidomain_design_namelist=./multidomain_design_control/multidomain_kalbarri2coralbay_historicaltsunamis.nml

# Model name
model_name=extremeTestCase_kalbarri2coralbay_B_tidaladjustment_DEBUG_theta19highres
# Ambient sea level of 0.0m AHD, approximation of expected sea level
ambient_sea_level=0

# Run the model
mpiexec -n 16 -map-by numa:SPAN -bind-to numa -x OMP_NUM_THREADS=13 -x OMP_PROC_BIND=true \
    ./model_debug \
    $stage_file \
    $multidomain_design_namelist \
    $model_name \
    full \
    $ambient_sea_level > outfile.log
