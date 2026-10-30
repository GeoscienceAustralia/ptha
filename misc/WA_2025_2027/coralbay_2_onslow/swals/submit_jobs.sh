#!bin/bash

for file in ./run_failed_job_low_ts_*.pbs; do
  qsub $file
  mv $file ./log/
done
  
