#!/bin/bash

cd ..

# qsub all .pbs files in the current directory
jobs=$(ls run_ptha/*.pbs)
# don't sumbmit the template
jobs=$(echo $jobs | sed 's/template.pbs//g')


for job in $jobs
do
    echo $job
    qsub $job
    mv $job run_ptha/submitted_jobs/
done

cd run_ptha
