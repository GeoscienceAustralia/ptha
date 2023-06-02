#!/bin/bash
#PBS -P w85
#PBS -q normal
#PBS -l walltime=24:00:00
#PBS -lmem=190GB
#PBS -lncpus=48
#PBS -ljobfs=20GB
#PBS -l wd
#PBS -l storage=scratch/w85+gdata/w85

source R_400_NCI_modules.sh

# The SWALS output directory, underneath ../../swals/OUTPUTS, where we check rasters from
runname='ptha18-GreaterPerth-sealevel60cm-reviseddomain-highres/random_outerrisesunda/'
# Find the raster outputs from one run (we assume all model runs have the same number of domains).
testfile=$(ls '../../swals/OUTPUTS/'$runname'ptha18_random_scenarios'_*/raster*.tar | head -n1)
# Count the rasters in the multidomain outputs
maxrasts=$(tar --list -f $testfile | grep "depth_" | wc -w)

Rscript compute_exceedance_rates_for_threshold_stage_logic_tree_mean_newParallelPartition.R  $runname $maxrasts

# The SWALS output directory, underneath ../../swals/OUTPUTS, where we check rasters from
runname='ptha18-GreaterPerth-sealevel60cm-reviseddomain-highres/random_sunda2/'
# Find the raster outputs from one run (we assume all model runs have the same number of domains).
testfile=$(ls '../../swals/OUTPUTS/'$runname'ptha18_random_scenarios'_*/raster*.tar | head -n1)
# Count the rasters in the multidomain outputs
maxrasts=$(tar --list -f $testfile | grep "depth_" | wc -w)

Rscript compute_exceedance_rates_for_threshold_stage_logic_tree_mean_newParallelPartition.R  $runname $maxrasts
