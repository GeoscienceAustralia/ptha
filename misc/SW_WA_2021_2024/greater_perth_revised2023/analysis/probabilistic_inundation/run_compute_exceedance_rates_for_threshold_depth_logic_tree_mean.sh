#!/bin/bash
#PBS -P w85
#PBS -q normal
#PBS -l walltime=08:00:00
#PBS -lmem=190GB
#PBS -lncpus=48
#PBS -ljobfs=20GB
#PBS -l wd
#PBS -l storage=scratch/w85+gdata/w85

source R_431_NCI_modules.sh

#
# outerrisesunda calculations
#

# The SWALS output directory, underneath ../../swals/OUTPUTS, where we check rasters from
runname='../../swals/OUTPUTS/ptha18-GreaterPerth2023-sealevel60cm/random_outerrisesunda/'
# Find the raster outputs from one run (we assume all model runs have the same number of domains).
testfile=$(ls $runname'ptha18_random_scenarios'_*/raster*.tar | head -n1)
# Count the rasters in the multidomain outputs
maxrasts=$(tar --list -f $testfile | grep "depth_" | wc -w)

Rscript compute_exceedance_rates_at_logic_tree_mean.R  $runname $maxrasts depth 0.001

#
# sunda2 calculations
#

# The SWALS output directory, underneath ../../swals/OUTPUTS, where we check rasters from
runname='../../swals/OUTPUTS/ptha18-GreaterPerth2023-sealevel60cm/random_sunda2/'
# Find the raster outputs from one run (we assume all model runs have the same number of domains).
testfile=$(ls $runname'ptha18_random_scenarios'_*/raster*.tar | head -n1)
# Count the rasters in the multidomain outputs
maxrasts=$(tar --list -f $testfile | grep "depth_" | wc -w)

Rscript compute_exceedance_rates_at_logic_tree_mean.R  $runname $maxrasts depth 0.001
