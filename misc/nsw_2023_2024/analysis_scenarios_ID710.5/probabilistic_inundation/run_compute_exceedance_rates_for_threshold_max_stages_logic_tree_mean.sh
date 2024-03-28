#!/bin/bash
#PBS -P w85
#PBS -q normal
#PBS -l walltime=24:00:00
#PBS -lmem=190GB
#PBS -lncpus=48
#PBS -ljobfs=20GB
#PBS -l wd
#PBS -l storage=scratch/w85+gdata/w85

source R_431_NCI_modules.sh

# Array with SWALS output directory, underneath ../../swals/OUTPUTS, where we check rasters from
all_runnames=(\
"../../swals/OUTPUTS/ptha18-NSW2023b-ID710.5-sealevel110cm/random_alaskaaleutians/" \
"../../swals/OUTPUTS/ptha18-NSW2023b-ID710.5-sealevel110cm/random_kermadectonga2/" \
"../../swals/OUTPUTS/ptha18-NSW2023b-ID710.5-sealevel110cm/random_newhebrides2/" \
"../../swals/OUTPUTS/ptha18-NSW2023b-ID710.5-sealevel110cm/random_outerrise_kermadectonga/" \
"../../swals/OUTPUTS/ptha18-NSW2023b-ID710.5-sealevel110cm/random_outerrisenewhebrides/" \
"../../swals/OUTPUTS/ptha18-NSW2023b-ID710.5-sealevel110cm/random_outerrise_puysegur/" \
"../../swals/OUTPUTS/ptha18-NSW2023b-ID710.5-sealevel110cm/random_puysegur2/" \
"../../swals/OUTPUTS/ptha18-NSW2023b-ID710.5-sealevel110cm/random_solomon2/" \
"../../swals/OUTPUTS/ptha18-NSW2023b-ID710.5-sealevel110cm/random_southamerica/")

first_run=${all_runnames[0]}

# Find the raster outputs from one run (we assume all model runs have the same number of domains).
testfile=$(ls $first_run'ptha18_random_scenarios'_*/raster*.tar | head -n1)
# Count the rasters in the multidomain outputs
maxrasts=$(tar --list -f $testfile | grep "max_stage_" | wc -w)

# Do the calculations for each source
for runname in ${all_runnames[@]}; do
  echo 'RUNNAME:' $runname $maxrasts;
  Rscript compute_exceedance_rates_at_logic_tree_mean.R $runname $maxrasts max_stage 1.101 2.1 3.1 4.1 5.1 6.1 7.1 8.1 9.1 10.1;
done

