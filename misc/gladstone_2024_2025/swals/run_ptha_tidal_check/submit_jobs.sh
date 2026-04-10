#!/bin/bash
# Submit hazard runs pbs scripts to the cluster

cwd=$(pwd)
cd ..

# for each .pbs file in the current directory
for file in $cwd/*.pbs; do
    echo "To submit $file"
done

# confirm the submission
echo "Submitting all .pbs files in the current directory"
read -p "Press enter to continue"; echo

for file in $cwd/*.pbs; do
    echo "Submitted $file"
    qsub $file
    mv $file $cwd/submitted_jobs
done

cd $cwd
