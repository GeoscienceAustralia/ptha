#!/bin/bash
#PBS -P w85
#PBS -q normal
#PBS -l walltime=48:00:00
#PBS -lmem=64GB
#PBS -lncpus=16
#PBS -l wd

# List all subdirectories at the level of individual tsunami runs
all_subdirs=$(ls -d ./unit_source_tsunami/*)

source ./template/swals_modules_nci.sh

export OMP_NUM_THREADS=1 # Run each model in serial

counter=0
one=1

mybasedir=$(pwd)

# Loop over all subdirectories
for i in $all_subdirs; do
    # If checkFile exists, the job is already submitted, so we move on
    checkFile=$i'/running_flag.txt';
    if [ -f $checkFile ];
    then
	#echo $i
	#echo 'ALREADY ON'
        continue 
    else
	# The job has not been submitted, so we should submit it
        cd $i
        touch 'running_flag.txt'
	#echo $(pwd)
        #echo 'Could run this one';
        #source run_jagurs_raw.sh &
        #./java_example < model_namelist_java2.in > outfile.log &

	# See if we can bind to individual CPUs for better speed
        numactl -C $counter ./generic_model model_namelist.in > outfile.log &
        counter=$(($counter+$one));
        #echo $counter;
	cd $mybasedir
    fi

    # Once we have submitted 16 jobs, break from this loop.
    if [ $counter = 16 ];
    then
        break
    fi
done      

wait
