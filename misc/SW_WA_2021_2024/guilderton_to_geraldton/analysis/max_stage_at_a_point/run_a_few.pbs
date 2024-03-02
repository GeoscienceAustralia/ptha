#!/bin/bash
#PBS -P w85
## Use the copy queue because it has internet access
#PBS -q copyq
# With 1 CPU it takes about 40 mins for each job
#PBS -l walltime=01:00:00
#PBS -lmem=48GB
#PBS -lncpus=1
#PBS -l wd
#PBS -l storage=scratch/w85+gdata/w85+gdata/fj6
#PBS -J 5-14
#PBS -r y

source ../../swals/modules_R_431.sh

# APPROACH 1 - in parallel
# no point on copy queue since limited to 1 CPU
# module load nci-parallel/1.0.0a

# export ncores_per_task=4
# export ncores_per_numanode=12

# mpirun -np $((PBS_NCPUS/ncores_per_task)) --map-by ppr:$((ncores_per_numanode/ncores_per_task)):NUMA:PE=${ncores_per_task} nci-parallel --input-file run_a_few.in --timeout 4000

# APPROACH 2 - do in series
# while read command_line; do
#     # Do each command line-by-line. Cut down to the first 4 words
#     $command_line | cut -d ' ' -f 1-4
# done <run_a_few.in


# APPROACH 3 - use array job
# Only 4 finished in 4 hours using approach 2.
# So run the remaining 10 using an array job #PBS -J above
$(sed "${PBS_ARRAY_INDEX}q;d" run_a_few.in | cut -d ' ' -f 1-4)
