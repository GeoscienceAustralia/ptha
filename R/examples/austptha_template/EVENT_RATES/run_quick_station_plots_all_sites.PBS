#!/bin/bash
#PBS -P w85
#PBS -q normal 
#PBS -l walltime=24:00:00
#PBS -lmem=64GB
#PBS -lncpus=16
#PBS -l wd

source R_modules.sh

# Number of cores on a single raijin node
export MC_CORES=16

# The hazard points can be ordered by longitude.
# This controls the 'percentiles of longitude' that we run
# e.g. to run the lower 10%;
#    export lower_percentile=0
#    export upper_percentile=10
# e.g. to run the middle 10%;
#    export lower_percentile=45
#    export upper_percentile=55
export lower_percentile=0
export upper_percentile=10

# Start $MC_CORES workers, which will repeatedly execute plot jobs
for i in $(seq 1 $MC_CORES); do
    Rscript quick_station_plots_all_sites.R $lower_percentile $upper_percentile $i $MC_CORES > /dev/null &
done

