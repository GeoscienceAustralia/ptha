#!/bin/bash
#PBS -P w85
#PBS -q normalsr
#PBS -l walltime=0:30:00
#PBS -lmem=1000GB
#PBS -lncpus=208
#PBS -l wd
#PBS -l storage=scratch/w85+gdata/w85

source SWALS_ifx_modules_2025_llvm.sh

# Initial stage 
stage_file=../sources/test_sources/small_source/outerrisesunda_row_0000760_Mw_72_HS.tif

## With tidal adjustment
#multidomain_design_namelist=./multidomain_design_control/multidomain_kalbarri2coralbay_hazard.nml

## No tidal adjustment
multidomain_design_namelist=./multidomain_design_control/multidomain_kalbarri2coralbay_B_historicaltsunamis.nml

# Model name
model_name=smallTestCase_kalbarri2coralbay_B_notidaladjustment
# Ambient sea level of 0.0m AHD, approximation of expected sea level
ambient_sea_level=0

# Run the model
mpiexec -n 16 -map-by numa:SPAN -bind-to numa -x OMP_NUM_THREADS=13 -x OMP_PROC_BIND=true \
    ./model_build \
    $stage_file \
    $multidomain_design_namelist \
    $model_name \
    test_load_balance \
    $ambient_sea_level > outfile.log
