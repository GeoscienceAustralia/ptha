#!/bin/bash

# qsub all .pbs files in the current directory
jobs=$(ls *.pbs)
# don't sumbmit the template
jobs=$(echo $jobs | sed 's/template.pbs//g')

for job in $jobs
do
    echo $job
    qsub $job
    mv $job submitted_jobs/
done
