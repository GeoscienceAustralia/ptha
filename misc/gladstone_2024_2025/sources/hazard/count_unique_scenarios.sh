#!/bin/bash

dirs=ptha_batch/random_*/scenario_initial_conditions

for dir in $dirs; do
    ls $dir | wc -l
done
 